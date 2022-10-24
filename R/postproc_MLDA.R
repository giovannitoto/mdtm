# ---------------------------------------------------------------------------- #

#' Post-processing for Hashtag-LDA (CGS)
#'
#' @description
#' Compute the logarithm of the a posteriori distribution at each MCMC iteration and the posterior means of the parameters.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param result_folder A string specifying the folder in which the results of a CGS obtained using the function \code{\link{CGS_MicroblogLDA}} are saved.
#' @param postproc_file A string specifying the RDS file in which the results will be saved. Write the name of the file without \code{.RDS} at the end.
#' @param iterations A vector of integers between 1 and the length of the chain which specifies the iterations to consider for the computation of the posterior means. Default is all the iterations.
#' @param verbose Logical: if \code{TRUE}, print the value of the logarithm of the a posteriori distribution at each iteration. Default is \code{TRUE}.
#'
#' @seealso \code{\link{CGS_MicroblogLDA}}, \code{\link{pred_MicroblogLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
postproc_MicroblogLDA <- function(result_folder, postproc_file, iterations = NULL, verbose = TRUE) {
  # -------------------------------------------------------------------------- #
  if(dir.exists(result_folder)) {
    hyper <- readRDS(file.path(result_folder, "hyperparameters.RDS"))
  } else {
    stop(result_folder, "' does not exist.\n", sep="")
  }
  # -------------------------------------------------------------------------- #
  if(postproc_file == "hyperparameters") {
    stop("'hyperparameters.RDS' not valid: select another value for the 'postproc_file' argument.\n", sep="")
  }
  # -------------------------------------------------------------------------- #
  if(is.null(iterations)) iterations <- 1:hyper$iterations
  # -------------------------------------------------------------------------- #
  N_sum <- apply(hyper$N, 2, sum)
  Dusers <- sapply(1:hyper$U, function(uu) sum(hyper$doc_users==uu))
  b <- simplify2array(hyper$b)
  loglik_list <- rep(0, hyper$iterations)
  x_est <- rep(0, hyper$D)
  zstar_est <- matrix(0, nrow = hyper$D, ncol = length(iterations))
  lambda_est <- matrix(0, nrow = hyper$D, ncol = hyper$T)
  y_est <- list()
  for (k in 1:hyper$K) {
    y_est[[k]] <- matrix(0, nrow = hyper$D, ncol = max(hyper$N[, k]))
  }
  thetastar_est <- matrix(0, nrow = hyper$U, ncol = hyper$T)
  theta_est <- matrix(0, nrow = hyper$D, ncol = hyper$T)
  piT_est <- 0
  delta_est <- 0
  phi_est <- list(); phiB_est <- list(); pi_est <- list()
  for (k in 1:hyper$K) {
    phi_est[[k]] <- matrix(0, nrow = hyper$T, ncol = hyper$V[[k]])
    phiB_est[[k]] <- rep(0, hyper$V[[k]])
    pi_est[[k]] <- 0
  }
  # -------------------------------------------------------------------------- #
  if(verbose) cat(as.character(Sys.time()), " Log-likelihood:\n", sep="")
  for (m in iterations) {
    # import m-th state of the chain
    x <- readRDS(file.path(result_folder, m, "x.RDS"))
    zstar <- readRDS(file.path(result_folder, m, "zstar.RDS"))
    lambda <- readRDS(file.path(result_folder, m, "lambda.RDS"))
    y <- list(); z <- list()
    for (k in 1:hyper$K) {
      y[[k]] <- readRDS(file.path(result_folder, m, paste("y", k, ".RDS", sep="")))
      z[[k]] <- readRDS(file.path(result_folder, m, paste("z", k, ".RDS", sep="")))
      if(m %in% iterations) y_est[[k]] <- y_est[[k]] + y[[k]]
    }
    if(m %in% iterations) x_est <- x_est + x
    if(m %in% iterations) zstar_est[, m] <- zstar
    if(m %in% iterations) lambda_est <- lambda_est + lambda
    # create matrices of counts
    X1 <- rep(0, hyper$U)
    Zstar <- matrix(0, nrow = hyper$U, ncol = hyper$T);
    LAMBDA1 = 0
    Z <- matrix(0, nrow = hyper$D, ncol = hyper$T);
    Yv1 <- rep(0, hyper$K);
    WY1ZX <- list(); WY0 <- list()
    for (k in 1:hyper$K) {
      WY1ZX[[k]] = matrix(0, nrow = hyper$V[k], ncol = hyper$T)
      WY0[[k]] = rep(0, hyper$V[k])
    }
    # update counts
    update_counts_MicroblogLDA(hyper$w, hyper$doc_users-1, Dusers, hyper$alphastar,
                               hyper$alpha, hyper$beta, b, hyper$bdelta,
                               hyper$bT, hyper$alpha0, hyper$T, hyper$K, hyper$U,
                               hyper$D, hyper$N, x, zstar, lambda, y, z,
                               X1, Zstar, LAMBDA1, Z, Yv1, WY1ZX, WY0, FALSE);
    # estimate thetastar
    thetastar <- hyper$alphastar + t(Zstar); thetastar <- t(thetastar)
    thetastar <- thetastar / apply(thetastar, 1, sum)
    if(m %in% iterations) thetastar_est <- thetastar_est + thetastar
    # estimate theta
    theta <- hyper$alpha + t(Z); theta <- t(theta)
    theta <- theta / apply(theta, 1, sum)
    if(m %in% iterations) theta_est <- theta_est + theta
    # estimate piT
    piT <- (hyper$bT[1] + X1) / sum(hyper$bT, Dusers)
    if(m %in% iterations) piT_est <- piT_est + piT
    # estimate delta
    delta <- (hyper$bdelta[1] + LAMBDA1)/ sum(hyper$bdelta, hyper$D * hyper$T)
    if(m %in% iterations) delta_est <- delta_est + delta
    # estimate phi, phiB, pi
    phi <- list(); phiB <- list(); Pi <- list()
    for (k in 1:hyper$K) {
      # estimate phi
      phi[[k]] <- hyper$beta[[k]] + WY1ZX[[k]]; phi[[k]] <- t(phi[[k]])
      phi[[k]] <- phi[[k]] / apply(phi[[k]], 1, sum)
      if(m %in% iterations) phi_est[[k]] <- phi_est[[k]] + phi[[k]]
      # estimate phiB
      phiB[[k]] <- (hyper$beta[[k]] + WY0[[k]]); phiB[[k]] <- phiB[[k]] / sum(phiB[[k]])
      if(m %in% iterations) phiB_est[[k]] <- phiB_est[[k]] + phiB[[k]]
      # estimate pi
      Pi[[k]] <- (hyper$b[[k]][1] + Yv1[[k]]) / sum(hyper$b[[k]], N_sum[k])
      if(m %in% iterations) pi_est[[k]] <- pi_est[[k]] + Pi[[k]]
    }
    # compute log-likelihood
    loglik_list[m] <- loglik_MicroblogLDA(hyper$alphastar, hyper$alpha,
                              hyper$alpha0, hyper$beta, hyper$bdelta, hyper$b, hyper$bT,
                              thetastar, theta, phi, phiB, delta, Pi, piT,
                              x, zstar, lambda, y, z,
                              X1, Zstar, LAMBDA1, Z, Yv1, WY1ZX, WY0,
                              hyper$D, hyper$T, Dusers, N_sum, hyper$K)
    if(verbose) cat(as.character(Sys.time()), "  - iteration ", m, ": ", loglik_list[m], "\n", sep="")
  }
  # MCMC estimates
  x_est <- x_est / length(iterations)
  zstar_est <- apply(zstar_est, 1, Mode)
  lambda_est <-lambda_est / length(iterations)
  for (k in 1:hyper$K) {
    y_est[[k]] <- y_est[[k]] / length(iterations)
  }
  thetastar_est <- thetastar_est / length(iterations)
  theta_est <- theta_est / length(iterations)
  lambda_est <- lambda_est / length(iterations)
  piT_est <- piT_est / length(iterations)
  delta_est <- delta_est / length(iterations)
  for (k in 1:hyper$K) {
    phi_est[[k]] <- phi_est[[k]] / length(iterations)
    phiB_est[[k]] <- phiB_est[[k]] / length(iterations)
    pi_est[[k]] <- pi_est[[k]] / length(iterations)
  }
  # -------------------------------------------------------------------------- #
  hyper[["loglik"]] <- loglik_list
  hyper[["iterations"]] <- iterations
  hyper[["x"]] <- x_est
  hyper[["zstar"]] <- zstar_est
  hyper[["lambda"]] <- lambda_est
  hyper[["y"]] <- y_est
  hyper[["piT"]] <- piT_est
  hyper[["thetastar"]] <- thetastar_est
  hyper[["theta"]] <- theta_est
  hyper[["delta"]] <- delta_est
  hyper[["phi"]] <- phi_est
  hyper[["phiB"]] <- phiB_est
  hyper[["pi"]] <- pi_est
  # -------------------------------------------------------------------------- #
  saveRDS(hyper, file.path(result_folder, paste(postproc_file, ".RDS", sep="")))
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
