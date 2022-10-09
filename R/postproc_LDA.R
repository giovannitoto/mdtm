# ---------------------------------------------------------------------------- #

postproc_LDA <- function(result_folder, postproc_file, iterations = NULL, verbose = TRUE) {
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
  loglik_list <- rep(0, hyper$iterations)
  theta_est <- matrix(0, nrow = hyper$D, ncol = hyper$T)
  phi_est <- matrix(0, nrow = hyper$T, ncol = hyper$V)
  # -------------------------------------------------------------------------- #
  if(verbose) cat(as.character(Sys.time()), " Log-likelihood:\n", sep="")
  for (m in iterations) {
    # import m-th state of the chain
    zV <- readRDS(file.path(result_folder, m, "zV.RDS"))
    # create matrices of counts
    WY1ZX <- matrix(0, nrow = hyper$V, ncol = hyper$T)
    Z <- matrix(0, nrow = hyper$D, ncol = hyper$T)
    # update counts
    update_counts_LDA(hyper$w, hyper$alpha, hyper$T, hyper$D, hyper$N,
                      zV, WY1ZX, Z, FALSE);
    # estimate theta
    theta <- hyper$alpha + t(Z); theta <- t(theta)
    theta <- theta / apply(theta, 1, sum)
    if(m %in% iterations) theta_est <- theta_est + theta
    # estimate phi
    phi <- hyper$betaV + WY1ZX; phi <- t(phi)
    phi <- phi / apply(phi, 1, sum)
    if(m %in% iterations) phi_est <- phi_est + phi
    # compute log-likelihood
    loglik_list[m] <- loglik_LDA(hyper$alpha, hyper$betaV, theta, phi, zV, WY1ZX, Z)
    if(verbose) cat(as.character(Sys.time()), "  - iteration ", m, ": ", loglik_list[m], "\n", sep="")
  }
  # MCMC estimates
  theta_est <- theta_est / length(iterations)
  phi_est <- phi_est / length(iterations)
  # -------------------------------------------------------------------------- #
  hyper[["loglik"]] <- loglik_list
  hyper[["iterations"]] <- iterations
  hyper[["theta"]] <- theta_est
  hyper[["phi"]] <- phi_est
  # -------------------------------------------------------------------------- #
  saveRDS(hyper, file.path(result_folder, paste(postproc_file, ".RDS", sep="")))
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
