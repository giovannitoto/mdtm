# ---------------------------------------------------------------------------- #

pred_LDA <- function(w, betaV_new = NULL, postproc_file, single_doc = TRUE,
                     iterations = 300, seed = 28, result_folder) {
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
  hyper <- list("w" = w, "alpha" = postproc$alpha, "betaV" = postproc$betaV,
                "iterations" = iterations, "seed" = seed, "T" = postproc$T,
                "D" = D, "V" = postproc$V, "N" = N, "single_doc" = single_doc,
                "postproc" = postproc_file)
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
  saveRDS(hyper, file.path(result_folder, "hyperparameters.RDS"))
  # -------------------------------------------------------------------------- #
  # get matrices of counts from the parameters
  Z <- postproc$theta * (sum(hyper$alpha) + sum(hyper$N))
  Z <- t(t(Z) - hyper$alpha)
  WY1ZX <- postproc$phi * (sum(hyper$betaV) + apply(Z, 2, sum))
  WY1ZX <- t(WY1ZX) - hyper$betaV
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder", "WY1ZX")))
  # -------------------------------------------------------------------------- #
  if(hyper$single_doc) {
    # update one new document at a time
    pred_single_LDA(hyper$w, hyper$alpha, hyper$betaV, hyper$iterations, hyper$T,
                    hyper$D, hyper$V, hyper$N, WY1ZX, result_folder)
  } else {
    # update all the new documents at the same time
    pred_all_LDA(hyper$w, hyper$alpha, hyper$betaV, hyper$iterations, hyper$T,
                 hyper$D, hyper$V, hyper$N, WY1ZX, result_folder)
  }
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}
