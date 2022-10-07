
postproc_HashtagLDA <- function(result_folder, postproc_file, iterations = NULL, verbose = TRUE) {
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
  L_sum <- sum(hyper$L)
  loglik_list <- rep(0, hyper$iterations)
  zstar_est <- matrix(0, nrow = hyper$D, ncol = length(iterations))
  thetastar_est <- matrix(0, nrow = hyper$U, ncol = hyper$T)
  phi_est <- matrix(0, nrow = hyper$T, ncol = hyper$V)
  psi_est <- matrix(0, nrow = hyper$T, ncol = hyper$H)
  psiB_est <- rep(0, hyper$H)
  piH_est <- 0
  # -------------------------------------------------------------------------- #
  if(verbose) cat(as.character(Sys.time()), " Log-likelihood:\n", sep="")
  for (m in iterations) {
    # import m-th state of the chain
    zstar <- readRDS(file.path(result_folder, m, "zstar.RDS"))
    yH <- readRDS(file.path(result_folder, m, "yH.RDS"))
    if(m %in% iterations) zstar_est[, m] <- zstar
    # create matrices of counts
    WY1ZX <- matrix(0, nrow = hyper$V, ncol = hyper$T)
    HY1ZX <- matrix(0, nrow = hyper$H, ncol = hyper$T)
    Zstar <- matrix(0, nrow = hyper$U, ncol = hyper$T)
    Yh1 <- 0
    HY0 <- rep(0, hyper$H)
    # update counts
    update_counts_HashtagLDA(hyper$w, hyper$h, hyper$doc_users-1, hyper$alphastar,
                             hyper$bH, hyper$T, hyper$D, hyper$N, hyper$L,
                             zstar, yH, WY1ZX, HY1ZX, Zstar, Yh1, HY0, FALSE);
    # estimate thetastar
    thetastar <- hyper$alphastar + t(Zstar); thetastar <- t(thetastar)
    thetastar <- thetastar / apply(thetastar, 1, sum)
    if(m %in% iterations) thetastar_est <- thetastar_est + thetastar
    # estimate phi
    phi <- hyper$betaV + WY1ZX; phi <- t(phi)
    phi <- phi / apply(phi, 1, sum)
    if(m %in% iterations) phi_est <- phi_est + phi
    # estimate psi
    psi <- hyper$betaH + HY1ZX; psi <- t(psi)
    psi <- psi / apply(psi, 1, sum)
    if(m %in% iterations) psi_est <- psi_est + psi
    # estimate psiB
    psiB <- (hyper$betaH + HY0); psiB <- psiB / sum(psiB)
    if(m %in% iterations) psiB_est <- psiB_est + psiB
    # estimate piH
    piH <- (hyper$bH[1] + Yh1) / sum(hyper$bH, L_sum)
    if(m %in% iterations) piH_est <- piH_est + piH
    # compute log-likelihood
    loglik_list[m] <- loglik_HashtagLDA(hyper$alphastar, hyper$betaV, hyper$betaH,
                                        hyper$bH, thetastar, phi, psi, psiB, piH,
                                        zstar, yH, WY1ZX, HY1ZX, Zstar, Yh1, HY0, L_sum)
    if(verbose) cat(as.character(Sys.time()), "  - iteration ", m, ": ", loglik_list[m], "\n", sep="")
  }
  # MCMC estimates
  zstar_est <- apply(zstar_est, 1, Mode)
  thetastar_est <- thetastar_est / length(iterations)
  phi_est <- phi_est / length(iterations)
  psi_est <- psi_est / length(iterations)
  psiB_est <- psiB_est / length(iterations)
  piH_est <- piH_est / length(iterations)
  # -------------------------------------------------------------------------- #
  output <- list("loglik" = loglik_list, "iterations" = iterations,
                 "thetastar" = thetastar_est, "phi" = phi_est, "psi" = psi_est,
                 "psiB" = psiB_est, "piH" = piH_est)
  saveRDS(output, file.path(result_folder, paste(postproc_file, ".RDS", sep="")))
  # -------------------------------------------------------------------------- #
}
