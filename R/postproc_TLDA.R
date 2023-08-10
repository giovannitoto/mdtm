# ---------------------------------------------------------------------------- #

#' Post-processing for Twitter-LDA (CGS)
#'
#' @description
#' Compute the logarithm of the a posteriori distribution at each MCMC iteration and the posterior means of the parameters.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param result_folder A string specifying the folder in which the results of a CGS obtained using the function \code{\link{CGS_TwitterLDA}} are saved.
#' @param postproc_file A string specifying the RDS file in which the results will be saved. Write the name of the file without \code{.RDS} at the end.
#' @param iterations A vector of integers between 1 and the length of the chain which specifies the iterations to consider for the computation of the posterior means. Default is all the iterations.
#' @param verbose Logical: if \code{TRUE}, print the value of the logarithm of the a posteriori distribution at each iteration. Default is \code{TRUE}.
#'
#' @seealso \code{\link{CGS_TwitterLDA}}, \code{\link{pred_TwitterLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
postproc_TwitterLDA <- function(result_folder, postproc_file, iterations = NULL, verbose = TRUE) {
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
  N_sum <- sum(hyper$N)
  loglik_list <- rep(0, hyper$iterations)
  zstar_est <- matrix(0, nrow = hyper$D, ncol = hyper$T)
  yV_est <- matrix(0, nrow = hyper$D, ncol = max(hyper$N))
  thetastar_est <- matrix(0, nrow = hyper$U, ncol = hyper$T)
  phi_est <- matrix(0, nrow = hyper$T, ncol = hyper$V)
  phiB_est <- rep(0, hyper$V)
  piV_est <- 0
  # -------------------------------------------------------------------------- #
  if(verbose) cat(as.character(Sys.time()), " Log-likelihood:\n", sep="")
  for (m in hyper$iterations) {
    # import m-th state of the chain
    zstar <- readRDS(file.path(result_folder, m, "zstar.RDS"))
    yV <- readRDS(file.path(result_folder, m, "yV.RDS"))
    if(m %in% iterations) zstar_est[cbind(1:hyper$D, zstar)] <- zstar_est[cbind(1:hyper$D, zstar)] + 1
    if(m %in% iterations) yV_est <- yV_est + yV
    # create matrices of counts
    WY1ZX <- matrix(0, nrow = hyper$V, ncol = hyper$T)
    Zstar <- matrix(0, nrow = hyper$U, ncol = hyper$T)
    Yv1 <- 0
    WY0 <- rep(0, hyper$V)
    # update counts
    update_counts_TwitterLDA(hyper$w, hyper$doc_users-1, hyper$alphastar,
                             hyper$bV, hyper$T, hyper$D, hyper$N,
                             zstar, yV, WY1ZX, Zstar, Yv1, WY0, FALSE);
    # estimate thetastar
    thetastar <- hyper$alphastar + t(Zstar); thetastar <- t(thetastar)
    thetastar <- thetastar / apply(thetastar, 1, sum)
    if(m %in% iterations) thetastar_est <- thetastar_est + thetastar
    # estimate phi
    phi <- hyper$betaV + WY1ZX; phi <- t(phi)
    phi <- phi / apply(phi, 1, sum)
    if(m %in% iterations) phi_est <- phi_est + phi
    # estimate phiB
    phiB <- (hyper$betaV + WY0); phiB <- phiB / sum(phiB)
    if(m %in% iterations) phiB_est <- phiB_est + phiB
    # estimate piV
    piV <- (hyper$bV[1] + Yv1) / sum(hyper$bV, N_sum)
    if(m %in% iterations) piV_est <- piV_est + piV
    # compute log-likelihood
    loglik_list[m] <- loglik_TwitterLDA(hyper$alphastar, hyper$betaV, hyper$bV,
                                        thetastar, phi, phiB, piV,
                                        zstar, yV, Zstar, WY1ZX, Yv1, WY0, N_sum)
    if(verbose) cat(as.character(Sys.time()), "  - iteration ", m, ":\t", loglik_list[m], "\n", sep="")
  }
  # MCMC estimates
  zstar_est <- apply(zstar_est, 1, which.max)
  yV_est <- yV_est / length(iterations)
  thetastar_est <- thetastar_est / length(iterations)
  phi_est <- phi_est / length(iterations)
  phiB_est <- phiB_est / length(iterations)
  piV_est <- piV_est / length(iterations)
  # -------------------------------------------------------------------------- #
  hyper[["loglik"]] <- loglik_list
  hyper[["iterations"]] <- iterations
  hyper[["zstar"]] <- zstar_est
  hyper[["yV"]] <- yV_est
  hyper[["thetastar"]] <- thetastar_est
  hyper[["phi"]] <- phi_est
  hyper[["phiB"]] <- phiB_est
  hyper[["piV"]] <- piV_est
  # -------------------------------------------------------------------------- #
  saveRDS(hyper, file.path(result_folder, paste(postproc_file, ".RDS", sep="")))
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
