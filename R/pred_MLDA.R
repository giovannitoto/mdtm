# ---------------------------------------------------------------------------- #

#' Prediction of new documents using Microblog-LDA (CGS)
#'
#' @description
#' Predict new documents using an already-estimated Microblog-LDA.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param w A \eqn{D\times N_{\text{max}}} matrix \eqn{\mathbf{w}} of integers in \eqn{\{1,\ldots,V+V_{\text{new}}\}}.
#' @param h A \eqn{D\times L_{\text{max}}} matrix \eqn{\mathbf{h}} of integers in \eqn{\{1,\ldots,H+H_{\text{new}}\}}.
#' @param doc_users A \eqn{D}-dimensional vector of integers specifying the single author of each document. The number of authors is set to \code{U = max(doc_users)}.
#' @param betaV_new A \eqn{V_{\text{new}}}-dimensional vector \eqn{\bm{\beta}_{\text{new}}^V} of positive numbers.
#' @param betaH_new A \eqn{H_{\text{new}}}-dimensional vector \eqn{\bm{\beta}_{\text{new}}^H} of positive numbers.
#' @param postproc_file A string specifying the RDS file in which the results of post-processing obtained using the function \code{\link{postproc_MicroblogLDA}} are saved.
#' @param single_doc Logical: if \code{TRUE}, the algorithm predicts one document at a time, otherwise all at once. Default is \code{TRUE}.
#' @param iterations An integer number of iterations. Default is 300.
#' @param seed Seed. Default is 28.
#' @param result_folder A string specifying the folder in which the results will be saved.
#'
#' @seealso \code{\link{CGS_MicroblogLDA}}, \code{\link{pred_MicroblogLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
pred_MicroblogLDA <- function(w, doc_users, beta_new = NULL, postproc_file,
                              single_doc = TRUE, iterations = 300,
                              seed = 28, result_folder) {
  # -------------------------------------------------------------------------- #

  # CHECK INPUTS HERE

  # -------------------------------------------------------------------------- #
  if(file.exists(postproc_file)) {
    postproc <- readRDS(postproc_file)
  } else {
    stop("'", postproc_file, "' does not exist.\n", sep="")
  }
  # -------------------------------------------------------------------------- #
  # Creo cartella in cui salvare gli stati della catena
  result_folder <- file.path(getwd(), result_folder)
  if(!dir.exists(result_folder)) {
    dir.create(result_folder)
    for (m in 0:iterations) {
      dir.create(file.path(result_folder, m))
    }
  } else {
    stop("'", result_folder, "' already exists: select another value for the 'result_folder' argument.\n", sep="")
  }
  # -------------------------------------------------------------------------- #
  # Fisso il seme
  set.seed(seed)
  # Definisco alcune quantita' utili
  D <- length(doc_users)
  N <- matrix(0, nrow = D, ncol = postproc$K)
  for (k in 1:K) {
    N[, k] <- apply(w[[k]], 1, function(x) sum(x > 0))
  }
  # -------------------------------------------------------------------------- #
  hyper <- list("w" = w, "doc_users" = doc_users, "alphastar" = postproc$alphastar,
                "alpha" = postproc$alpha, "beta" = postproc$beta, "b" = postproc$b,
                "bdelta" = postproc$bdelta, "bT" = postproc$bT, "alpha0" = postproc$alpha0,
                "iterations" = iterations, "seed" = seed, "T" = postproc$T,
                "K" = postproc$K, "U" = postproc$U, "D" = D, "V" = postproc$V,
                "N" = N, "single_doc" = single_doc, "postproc" = postproc_file)
  rm(list = setdiff(ls(), c("beta_new", "hyper", "postproc", "result_folder")))
  for (k in 1:hyper$K) {
    new_words <- max(hyper$w[[k]]) - hyper$V[k]
    if(new_words > 0) {
      if(is.null(beta_new[[k]]) || length(beta_new[[k]]) != new_words) {
        stop(paste("'beta_new[[", k, "]]' not valid: must be a ", new_words, "-dimensional vector of positive numbers.", sep=""))
      } else {
        # update size of the vocabulary
        hyper$beta[[k]] <- c(hyper$beta[[k]], betaV_new[[k]])
        hyper$V[k] <- hyper$V[k] + new_words
        # add empty columns to phi
        postproc$phi[[k]] <- cbind(postproc$phi[[k]], matrix(0, nrow = hyper$T, ncol = new_words))
      }
    }
  }
  new_users <- max(hyper$doc_users) - hyper$U
  if(new_users > 0) {
    # update size of the vocabulary
    hyper$U <- hyper$U + new_users
    # add empty rows to thetastar
    postproc$thetastar <- cbind(postproc$thetastar, matrix(hyper$alphastar, nrow = new_users, ncol = hyper$T, byrow = TRUE))
  }
  saveRDS(hyper, file.path(result_folder, "hyperparameters.RDS"))
  # -------------------------------------------------------------------------- #
  # get matrices of counts from the posterior estimates
  X1 <- rep(0, hyper$U)
  Zstar <- matrix(0, nrow = hyper$U, ncol = hyper$T);
  LAMBDA1 <- 0
  Z <- matrix(0, nrow = postproc$D, ncol = hyper$T);
  Yv1 <- rep(0, hyper$K);
  WY1ZX <- list(); WY0 <- list()
  for (k in 1:hyper$K) {
    WY1ZX[[k]] = matrix(0, nrow = hyper$V[k], ncol = hyper$T)
    WY0[[k]] = rep(0, hyper$V[k])
  }
  # update counts
  b <- simplify2array(hyper$b)
  update_counts_MicroblogLDA(postproc$w, postproc$doc_users-1, postproc$alphastar,
                             postproc$alpha, postproc$beta, b, postproc$bdelta,
                             postproc$bT, postproc$alpha0, postproc$T, postproc$K,
                             postproc$D, postproc$N, postproc$x, postproc$zstar,
                             postproc$lambda, postproc$y, postproc$z,
                             X1, Zstar, LAMBDA1, Z, Yv1, WY1ZX, WY0, FALSE);
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "b", "result_folder", "WY1ZX", "HY1ZX", "Zstar", "Yh1", "HY0")))
  # -------------------------------------------------------------------------- #
  N_sum <- apply(postproc$N, 2, sum) + apply(hyper$N, 2, sum)
  Dusers <- sapply(1:hyper$U, function(uu) sum(postproc$doc_users == uu))
  Dusers <- Dusers + sapply(1:hyper$U, function(uu) sum(hyper$doc_users == uu))
  if(hyper$single_doc) {
    # update one new document at a time
    pred_single_MicroblogLDA(hyper$w, hyper$doc_users-1, hyper$alphastar, hyper$alpha,
                             hyper$beta, b, hyper$bdelta, hyper$bT, hyper$alpha0,
                             hyper$iterations, hyper$T, hyper$K, hyper$U, hyper$D,
                             hyper$V, hyper$N, N_sum, Dusers, X1, Zstar, LAMBDA1,
                             Yv1, WY1ZX, WY0, result_folder)
  } else {
    # update all the new documents at the same time
    pred_all_MicroblogLDA(hyper$w, hyper$doc_users-1, hyper$alphastar, hyper$alpha,
                             hyper$beta, b, hyper$bdelta, hyper$bT, hyper$alpha0,
                             hyper$iterations, hyper$T, hyper$K, hyper$U, hyper$D,
                             hyper$V, hyper$N, N_sum, Dusers, X1, Zstar, LAMBDA1,
                             Yv1, WY1ZX, WY0, result_folder)
  }
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}
