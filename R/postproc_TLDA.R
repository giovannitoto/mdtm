
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
  zstar_est <- matrix(0, nrow = hyper$D, ncol = length(iterations))
  thetastar_est <- matrix(0, nrow = hyper$U, ncol = hyper$T)
  phi_est <- matrix(0, nrow = hyper$T, ncol = hyper$V)
  phiB_est <- rep(0, hyper$V)
  piV_est <- 0
  # -------------------------------------------------------------------------- #
  if(verbose) cat(as.character(Sys.time()), " Log-likelihood:\n", sep="")
  for (m in iterations) {
    # import m-th state of the chain
    zstar <- readRDS(file.path(result_folder, m, "zstar.RDS"))
    yV <- readRDS(file.path(result_folder, m, "yV.RDS"))
    if(m %in% iterations) zstar_est[, m] <- zstar
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
    if(verbose) cat(as.character(Sys.time()), "  - iteration ", m, ": ", loglik_list[m], "\n", sep="")
  }
  # MCMC estimates
  zstar_est <- apply(zstar_est, 1, Mode)
  thetastar_est <- thetastar_est / length(iterations)
  phi_est <- phi_est / length(iterations)
  phiB_est <- phiB_est / length(iterations)
  piV_est <- piV_est / length(iterations)
  # -------------------------------------------------------------------------- #
  output <- list("loglik" = loglik_list, "iterations" = iterations,
                 "thetastar" = thetastar_est, "phi" = phi_est, "phiB" = phiB_est,
                 "piV" = piV)
  saveRDS(output, file.path(result_folder, paste(postproc_file, ".RDS", sep="")))
  # -------------------------------------------------------------------------- #
}
