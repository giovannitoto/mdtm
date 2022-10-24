# ---------------------------------------------------------------------------- #

#' Hashtag-LDA (CGS)
#'
#' @description
#' Implementation of the Collapsed Gibbs Sampler (CGS) for Hashtag-LDA.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param w A \eqn{D\times N_{\text{max}}} matrix \eqn{\mathbf{w}} of integers in \eqn{\{1,\ldots,V\}}.
#' @param h A \eqn{D\times L_{\text{max}}} matrix \eqn{\mathbf{h}} of integers in \eqn{\{1,\ldots,H\}}.
#' @param doc_users A \eqn{D}-dimensional vector of integers specifying the single author of each document. The number of authors is set to \code{U = max(doc_users)}.
#' @param alphastar A \eqn{T}-dimensional vector of positive numbers, \eqn{\bm{\alpha}^*}. The length of this vector specifies the number of topics \eqn{T}.
#' @param betaV A \eqn{V}-dimensional vector \eqn{\bm{\beta}^V} of positive numbers.
#' @param betaH A \eqn{H}-dimensional vector \eqn{\bm{\beta}^H} of positive numbers.
#' @param bH A \eqn{2}-dimensional vector \eqn{\bm{b}^H} of positive numbers.
#' @param iterations An integer number of iterations. Default is 300.
#' @param seed Seed. Default is 28.
#' @param result_folder A string specifying the folder in which the results will be saved.
#'
#' @seealso \code{\link{postproc_HashtagLDA}}, \code{\link{pred_HashtagLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
CGS_HashtagLDA <- function(w, h, doc_users, alphastar, betaV, betaH, bH,
                           iterations=300, seed=28, result_folder) {
  # -------------------------------------------------------------------------- #
  # Argomenti della funzione:
  #          w : matrice D x Nmax   | n-ma parola del d-mo documento (1,...,V)
  #          h : matrice D x Lmax   | l-mo hashtag del d-mo documento (1,...,H)
  #  doc_users : vettore D x 1      | autore del d-mo topic (1,...,U)
  #  alphastar : vettore TOPICS x 1 | parametro Dirichlet sul simplesso dei topic
  #      betaV : vettore V x 1      | parametro Dirichlet sul simplesso delle parole
  #      betaH : vettore H x 1      | parametro Dirichlet sul simplesso degli hashtag
  #         bH : vettore 2x1        | parametro Beta
  # iterations : intero             | numero di stati della catena da campionare
  #       seed : intero             | seme per rendere i risultati replicabili
  # -------------------------------------------------------------------------- #

  # CHECK INPUTS HERE

  # -------------------------------------------------------------------------- #
  cat(as.character(Sys.time()), " START\n\n", sep="")
  # Fisso il seme
  set.seed(seed)
  # Definisco alcune quantita' utili
  TOPICS <- length(alphastar)
  U <- max(doc_users)
  D <- length(doc_users)
  V <- length(betaV)
  H <- length(betaH)
  N <- apply(w, 1, function(x) sum(x > 0))
  L <- apply(h, 1, function(x) sum(x > 0))
  # -------------------------------------------------------------------------- #
  # Creo cartella in cui salvare gli stati della catena
  result_folder <- file.path(getwd(), result_folder)
  if(!dir.exists(result_folder)) {
    dir.create(result_folder)
    for (m in 0:iterations) {
      dir.create(file.path(result_folder, m))
    }
  } else {
    stop("\n'",result_folder,"' already exists: select another value for the 'result_folder' argument.\n", sep="")
  }
  # -------------------------------------------------------------------------- #
  hyper <- list("w" = w, "h" = h, "doc_users" = doc_users,
                "alphastar" = alphastar, "betaV" = betaV, "betaH" = betaH, "bH" = bH,
                "iterations" = iterations, "seed" = seed,
                "T" = TOPICS, "U" = U, "D" = D, "V" = V, "H" = H, "N" = N, "L" = L)
  saveRDS(hyper, file.path(result_folder, "hyperparameters.RDS"))
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder")))
  # -------------------------------------------------------------------------- #
  rcpp_CGS_HashtagLDA(hyper$w, hyper$h, hyper$doc_users-1, hyper$alphastar,
                      hyper$betaV, hyper$betaH, hyper$bH, hyper$iterations,
                      hyper$T, hyper$U, hyper$D, hyper$V, hyper$H, hyper$N,
                      hyper$L, result_folder)
  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " END", sep="")
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
