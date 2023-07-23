# ---------------------------------------------------------------------------- #

#' Latent Dirichlet Allocation (CGS)
#'
#' @description
#' Implementation of the Collapsed Gibbs Sampler (CGS) for Latent Dirichlet Allocarion.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param w A \eqn{D\times N_{\text{max}}} matrix \eqn{\mathbf{w}} of integers in \eqn{\{1,\ldots,V\}}.
#' @param alpha A \eqn{T}-dimensional vector \eqn{\bm{\alpha}} of positive numbers. The length of this vector specifies the number of topics \eqn{T}.
#' @param betaV A \eqn{V}-dimensional vector \eqn{\bm{\beta}^V} of positive numbers.
#' @param vocab A \eqn{V}-dimensional vector of vocabulary words.
#' @param iterations An integer number of iterations. Default is 300.
#' @param seed Seed. Default is 28.
#' @param result_folder A string specifying the folder in which the results will be saved.
#'
#' @seealso \code{\link{postproc_LDA}}, \code{\link{pred_LDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
CGS_LDA <- function(w, alpha, betaV, vocab = NULL, iterations = 300, seed = 28, result_folder) {
  # -------------------------------------------------------------------------- #
  # Argomenti della funzione:
  #          w : matrice D x Nmax   | n-ma parola del d-mo documento (1,...,V)
  #      alpha : vettore TOPICS x 1 | parametro Dirichlet sul simplesso dei topic
  #      betaV : vettore V x 1      | parametro Dirichlet sul simplesso delle parole
  # iterations : intero             | numero di stati della catena da campionare
  #       seed : intero             | seme per rendere i risultati replicabili
  # -------------------------------------------------------------------------- #
  if(!is.null(vocab)) {
    if(length(betaV) != length(vocab)) {
      stop("'vocab' not valid: it must be a ", length(betaV), "-dimensional vector.", sep="")
    }
  }

  # CHECK INPUTS HERE

  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " START\n\n", sep="")
  # Fisso il seme
  set.seed(seed)
  # Definisco alcune quantita' utili
  TOPICS <- length(alpha)
  D <- nrow(w)
  V <- length(betaV)
  N <- apply(w, 1, function(x) sum(x > 0))
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
  hyper <- list("w" = w, "alpha" = alpha, "betaV" = betaV,
                "iterations" = iterations, "seed" = seed,
                "T" = TOPICS, "D" = D, "V" = V, "N" = N, "vocab" = vocab)
  saveRDS(hyper, file.path(result_folder, "hyperparameters.RDS"))
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder")))
  # -------------------------------------------------------------------------- #
  rcpp_CGS_LDA(hyper$w, hyper$alpha, hyper$betaV, hyper$iterations, hyper$T,
               hyper$D, hyper$V, hyper$N, result_folder)
  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " END", sep="")
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
