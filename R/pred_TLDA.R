# ---------------------------------------------------------------------------- #

pred_TwitterLDA <- function(w, doc_users, betaV_new = NULL, postproc_file,
                            single_doc = TRUE, iterations = 300, seed = 28,
                            result_folder) {
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
  # -------------------------------------------------------------------------- #
  hyper <- list("w" = w, "doc_users" = doc_users, "alphastar" = postproc$alphastar,
                "betaV" = postproc$betaV, "bV" = postproc$bV, "iterations" = iterations,
                "seed" = seed, "T" = postproc$T, "U" = postproc$U, "D" = D, "V" = postproc$V,
                "N" = N, "single_doc" = single_doc, "postproc" = postproc_file)
  rm(list = setdiff(ls(), c("betaV_new", "hyper", "postproc", "result_folder")))
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
  ZY1 <- sapply(1:hyper$T, function(tt) sum(postproc$yV[postproc$zstar == tt, ]))
  WY1ZX <- postproc$phi * (sum(hyper$betaV) + ZY1)
  WY1ZX <- t(WY1ZX) - hyper$betaV
  Dusers <- sapply(1:hyper$U, function(uu) sum(hyper$doc_users == uu))
  Zstar <- postproc$theta * (sum(hyper$alphastar) + Dusers)
  Zstar <- t(t(Zstar) - hyper$alphastar)
  Yv1 <- postproc$piV * sum(hyper$bV, hyper$N) - hyper$bV[1]
  WY0 <- postproc$phiB * (sum(hyper$betaV, hyper$N) - Yv1) - hyper$betaV
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder", "WY1ZX", "Zstar", "Yv1", "WY0")))
  # -------------------------------------------------------------------------- #
  if(hyper$single_doc) {
    # update one new document at a time
    pred_single_TwitterLDA(hyper$w, hyper$doc_users-1, hyper$alphastar, hyper$betaV,
                           hyper$bV, hyper$iterations, hyper$T, hyper$U, hyper$D,
                           hyper$V, hyper$N, WY1ZX, Zstar, Yv1, WY0, result_folder)
  } else {
    # update all the new documents at the same time
    pred_all_TwitterLDA(hyper$w, hyper$doc_users-1, hyper$alphastar, hyper$betaV,
                        hyper$bV, hyper$iterations, hyper$T, hyper$U, hyper$D,
                        hyper$V, hyper$N, WY1ZX, Zstar, Yv1, WY0, result_folder)
  }
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}
