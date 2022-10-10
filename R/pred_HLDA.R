# ---------------------------------------------------------------------------- #

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
  D <- nrow(w)
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
  # get matrices of counts from the parameters
  ZY1 <- sapply(1:hyper$T, function(tt) sum(hyper$N[postproc$zstar == tt]))
  WY1ZX <- postproc$phi * (sum(hyper$betaV) + ZY1)
  WY1ZX <- t(WY1ZX) - hyper$betaV
  ZY1 <- sapply(1:hyper$T, function(tt) sum(postproc$yH[postproc$zstar == tt, ]))
  HY1ZX <- postproc$psi * (sum(hyper$betaH) + ZY1)
  HY1ZX <- t(HY1ZX) - hyper$betaH
  Dusers <- sapply(1:hyper$U, function(uu) sum(hyper$doc_users == uu))
  Zstar <- postproc$theta * (sum(hyper$alphastar) + Dusers)
  Zstar <- t(t(Zstar) - hyper$alphastar)
  Yh1 <- postproc$piH * sum(hyper$bH, hyper$L) - hyper$bH[1]
  HY0 <- postproc$psiB * (sum(hyper$betaH, hyper$L) - Yh1) - hyper$betaH
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder", "WY1ZX", "HY1ZX", "Zstar", "Yh1", "HY0")))
  # -------------------------------------------------------------------------- #
  if(hyper$single_doc) {
    # update one new document at a time
    pred_single_HashtagLDA(hyper$w, hyper$h, hyper$doc_users-1, hyper$alphastar,
                           hyper$betaV, hyper$betaH, hyper$bH, hyper$iterations,
                           hyper$T, hyper$U, hyper$D, hyper$V, hyper$H, hyper$N,
                           hyper$L, WY1ZX, HY1ZX, Zstar, Yh1, HY0, result_folder)
  } else {
    # update all the new documents at the same time
    pred_all_HashtagLDA(hyper$w, hyper$h, hyper$doc_users-1, hyper$alphastar,
                        hyper$betaV, hyper$betaH, hyper$bH, hyper$iterations,
                        hyper$T, hyper$U, hyper$D, hyper$V, hyper$H, hyper$N,
                        hyper$L, WY1ZX, HY1ZX, Zstar, Yh1, HY0, result_folder)
  }
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}
