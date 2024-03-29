# ---------------------------------------------------------------------------- #

#' Prediction of new documents using Hashtag-LDA (CGS)
#'
#' @description
#' Predict new documents using an already-estimated Hashtag-LDA.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param w A \eqn{D\times N_{\text{max}}} matrix \eqn{\mathbf{w}} of integers in \eqn{\{1,\ldots,V+V_{\text{new}}\}}.
#' @param h A \eqn{D\times L_{\text{max}}} matrix \eqn{\mathbf{h}} of integers in \eqn{\{1,\ldots,H+H_{\text{new}}\}}.
#' @param doc_users A \eqn{D}-dimensional vector of integers specifying the single author of each document. The number of authors is set to \code{U = max(doc_users)}.
#' @param betaV_new A \eqn{V_{\text{new}}}-dimensional vector \eqn{\bm{\beta}_{\text{new}}^V} of positive numbers.
#' @param betaH_new A \eqn{H_{\text{new}}}-dimensional vector \eqn{\bm{\beta}_{\text{new}}^H} of positive numbers.
#' @param postproc_file A string specifying the RDS file in which the results of post-processing obtained using the function \code{\link{postproc_HashtagLDA}} are saved.
#' @param single_doc Logical: if \code{TRUE}, the algorithm predicts one document at a time, otherwise all at once. Default is \code{TRUE}.
#' @param iterations An integer number of iterations. Default is 300.
#' @param seed Seed. Default is 28.
#' @param result_folder A string specifying the folder in which the results will be saved.
#'
#' @seealso \code{\link{CGS_HashtagLDA}}, \code{\link{pred_HashtagLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
pred_HashtagLDA <- function(w, h, doc_users, betaV_new = NULL, betaH_new = NULL,
                            postproc_file, single_doc = TRUE, iterations = 300,
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
  N <- apply(w, 1, function(x) sum(x > 0))
  L <- apply(h, 1, function(x) sum(x > 0))
  # -------------------------------------------------------------------------- #
  hyper <- list("w" = w, "h" = h, "doc_users" = doc_users, "alphastar" = postproc$alphastar,
                "betaV" = postproc$betaV, "betaH" = postproc$betaH, "bH" = postproc$bH,
                "iterations" = iterations, "seed" = seed, "T" = postproc$T, "U" = postproc$U,
                "D" = D, "V" = postproc$V, "H" = postproc$H, "N" = N, "L" = L,
                "single_doc" = single_doc, "postproc" = postproc_file)
  rm(list = setdiff(ls(), c("betaV_new", "betaH_new", "hyper", "postproc", "result_folder")))
  new_words <- max(hyper$w) - hyper$V
  if(new_words > 0) {
    if(is.null(betaV_new) || length(betaV_new) != new_words) {
      stop(paste("'betaV_new' not valid: must be a ", new_words, "-dimensional vector of positive numbers.", sep=""))
    } else {
      # update size of the vocabulary
      hyper$betaV <- c(hyper$betaV, betaV_new)
      hyper$V <- hyper$V + new_words
      # add empty columns to phi
      postproc$phi <- cbind(postproc$phi, matrix(0, nrow = hyper$T, ncol = new_words))
    }
  }
  new_hashtags <- max(hyper$h) - hyper$H
  if(new_hashtags > 0) {
    if(is.null(betaH_new) || length(betaH_new) != new_hashtags) {
      stop(paste("'betaH_new' not valid: must be a ", new_hashtags, "-dimensional vector of positive numbers.", sep=""))
    } else {
      # update size of the vocabulary
      hyper$betaH <- c(hyper$betaH, betaH_new)
      hyper$H <- hyper$H + new_hashtags
      # add empty columns to phi
      postproc$psi <- cbind(postproc$psi, matrix(0, nrow = hyper$T, ncol = new_hashtags))
      postproc$psiB <- c(postproc$psiB, rep(0, new_hashtags))
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
  WY1ZX <- matrix(0, nrow = hyper$V, ncol = hyper$T)
  HY1ZX <- matrix(0, nrow = hyper$H, ncol = hyper$T)
  Zstar <- matrix(0, nrow = hyper$U, ncol = hyper$T)
  Yh1 <- 0
  HY0 <- rep(0, hyper$H)
  # update counts
  update_counts_HashtagLDA(postproc$w, postproc$h, postproc$doc_users-1, postproc$alphastar,
                           postproc$bH, postproc$T, postproc$D, postproc$N, postproc$L,
                           postproc$zstar, postproc$yH, WY1ZX, HY1ZX, Zstar, Yh1, HY0, FALSE)
  # other quantities
  L_sum <- sum(postproc$L) + sum(hyper$L)
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "L_sum", "WY1ZX", "HY1ZX", "Zstar", "Yh1", "HY0", "result_folder")))
  # -------------------------------------------------------------------------- #
  if(hyper$single_doc) {
    # update one new document at a time
    pred_single_HashtagLDA(hyper$w, hyper$h, hyper$doc_users-1, hyper$alphastar,
                           hyper$betaV, hyper$betaH, hyper$bH, hyper$iterations,
                           hyper$T, hyper$U, hyper$D, hyper$V, hyper$H, hyper$N,
                           hyper$L, L_sum, WY1ZX, HY1ZX, Zstar, Yh1, HY0, result_folder)
  } else {
    # update all the new documents at the same time
    pred_all_HashtagLDA(hyper$w, hyper$h, hyper$doc_users-1, hyper$alphastar,
                        hyper$betaV, hyper$betaH, hyper$bH, hyper$iterations,
                        hyper$T, hyper$U, hyper$D, hyper$V, hyper$H, hyper$N,
                        hyper$L, L_sum, WY1ZX, HY1ZX, Zstar, Yh1, HY0, result_folder)
  }
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}
