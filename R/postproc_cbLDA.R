# ---------------------------------------------------------------------------- #

#' Post-processing for clickbait-LDA (CGS)
#'
#' @description
#' Compute the posterior means of the parameters.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param result_folder A string specifying the folder in which the results of a CGS obtained using the function \code{\link{CGS_TwitterLDA}} are saved.
#' @param postproc_file A string specifying the RDS file in which the results will be saved. Write the name of the file without \code{.RDS} at the end.
#' @param iterations A vector of integers between 1 and the length of the chain which specifies the iterations to consider for the computation of the posterior means. Default is all the iterations.
#' @param verbose Logical: if \code{TRUE}, print the value of the logarithm of the a posteriori distribution at each iteration. Default is \code{TRUE}.
#'
# @seealso \code{\link{CGS_TwitterLDA}}, \code{\link{pred_TwitterLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
postproc_clickbaitLDA <- function(result_folder, postproc_file, iterations = NULL, verbose = TRUE) {
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
  Ntot <- sum(hyper$N)
  x_est <- rep(0, hyper$D)
  zstar_est <- matrix(0, nrow = hyper$D, ncol = hyper$T)
  yV_est <- matrix(0, nrow = hyper$D, ncol = max(hyper$N))
  thetastar_est <- rep(0, hyper$T)
  phi_est <- matrix(0, nrow = hyper$T, ncol = hyper$V)
  phi_ncb_est <- rep(0, hyper$V)
  phi_cb_est <- rep(0, hyper$V)
  pi_ncb_est <- 0
  pi_cb_est <- 0
  pi_doc_est <- 0
  # -------------------------------------------------------------------------- #
  if(verbose) cat(as.character(Sys.time()), " Log-likelihood:\n", sep="")
  for (m in iterations) {
    # import m-th state of the chain
    x <- readRDS(file.path(result_folder, m, "x.RDS"))
    zstar <- readRDS(file.path(result_folder, m, "zstar.RDS"))
    yV <- readRDS(file.path(result_folder, m, "yV.RDS"))
    if(m %in% iterations) x_est <- x_est + x
    if(m %in% iterations) zstar_est[cbind(1:hyper$D, zstar)] <- zstar_est[cbind(1:hyper$D, zstar)] + 1
    if(m %in% iterations) yV_est <- yV_est + yV
    # create matrices of counts
    WY1ZX <- matrix(0, nrow = hyper$T, ncol = hyper$V)
    Zstar <- rep(0, hyper$T)
    X1 <- 0
    WY0 <- matrix(0, nrow = 2, ncol = hyper$V)
    Yv1 <- c(0,0)
    # update counts
    update_counts_clickbaitLDA(hyper$w, x, hyper$alphastar,
                               hyper$betaV, hyper$betaB, hyper$b_doc, hyper$b_back,
                               hyper$T, hyper$D, hyper$Dt, hyper$V, hyper$N,
                               zstar, yV, WY1ZX, Zstar, X1, Yv1, WY0, FALSE);
    # define counts
    N_sum0 <- sum(hyper$N * (1-x))
    N_sum1 <- sum(hyper$N * x)
    # estimate thetastar
    thetastar <- hyper$alphastar + Zstar
    thetastar <- thetastar / sum(thetastar)
    if(m %in% iterations) thetastar_est <- thetastar_est + thetastar
    # estimate phi
    phi <- hyper$betaV + WY1ZX
    phi <- phi / apply(phi, 1, sum)
    if(m %in% iterations) phi_est <- phi_est + phi
    # estimate phi_ncb
    phi_ncb <- hyper$betaB[2] + WY0[2, ]
    phi_ncb <- phi_ncb / sum(phi_ncb)
    if(m %in% iterations) phi_ncb_est <- phi_ncb_est + phi_ncb
    # estimate phi_cb
    phi_cb <- hyper$betaB[1] + WY0[1, ]
    phi_cb <- phi_cb / sum(phi_cb)
    if(m %in% iterations) phi_cb_est <- phi_cb_est + phi_cb
    # estimate pi_ncb
    pi_ncb <- (hyper$b_back[2,1] + Yv1[2]) / sum(hyper$b_back[2,], N_sum1)
    if(m %in% iterations) pi_ncb_est <- pi_ncb_est + pi_ncb
    # estimate pi_cb
    pi_cb <- (hyper$b_back[1,1] + Yv1[1]) / sum(hyper$b_back[1,], N_sum0)
    if(m %in% iterations) pi_cb_est <- pi_cb_est + pi_cb
    # estimate pi
    pi_doc <- (hyper$b_doc[1] + X1) / sum(hyper$b_doc, Ntot)
    if(m %in% iterations) pi_doc_est <- pi_doc_est + pi_doc
    # end of iteration
    if(verbose) cat(as.character(Sys.time()), "  - iteration ", m, "\n", sep="")
  }
  # MCMC estimates
  x_est <- x_est / length(iterations)
  zstar_est <- apply(zstar_est, 1, which.max)
  yV_est <- yV_est / length(iterations)
  thetastar_est <- thetastar_est / length(iterations)
  phi_est <- phi_est / length(iterations)
  phi_ncb_est <- phi_ncb_est / length(iterations)
  phi_cb_est <- phi_cb_est / length(iterations)
  pi_ncb_est <- pi_ncb_est / length(iterations)
  pi_cb_est <- pi_cb_est / length(iterations)
  pi_doc_est <- pi_doc_est / length(iterations)
  # -------------------------------------------------------------------------- #
  hyper[["iterations"]] <- iterations
  hyper[["x"]] <- x_est
  hyper[["zstar"]] <- zstar_est
  hyper[["yV"]] <- yV_est
  hyper[["thetastar"]] <- thetastar_est
  hyper[["phi"]] <- phi_est
  hyper[["phi_ncb"]] <- phi_ncb_est
  hyper[["phi_cb"]] <- phi_cb_est
  hyper[["pi_ncb"]] <- pi_ncb_est
  hyper[["pi_cb"]] <- pi_cb_est
  hyper[["pi_doc"]] <- pi_doc_est
  # -------------------------------------------------------------------------- #
  saveRDS(hyper, file.path(result_folder, paste(postproc_file, ".RDS", sep="")))
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
