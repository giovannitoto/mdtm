# ---------------------------------------------------------------------------- #

#' Click bait-LDA (CGS)
#'
#' @description
#' Implementation of the Collapsed Gibbs Sampler (CGS) for click bait-LDA (cb-LDA).
#'
#' @details
#' Implementation in R and C++.
#'
#' @param w A \eqn{D\times N_{\text{max}}} matrix \eqn{\mathbf{w}} of integers in \eqn{\{1,\ldots,V\}}.
#' @param x A \eqn{D^t}-dimensional vector \eqn{\bm{x}^t} of values in \eqn{\{0,1\}}, where \eqn{0leq D^t\leq D}.
#' @param alphastar A positive number \eqn{\alpha^*}.
#' @param betaV A positive number \eqn{\beta^V}.
#' @param beta_ncb A positive number \eqn{\beta^{ncb}}.
#' @param beta_cb A positive number \eqn{\beta^{cb}}.
#' @param bV A \eqn{2}-dimensional vector \eqn{\bm{b}^V} of positive numbers.
#' @param b_ncb A \eqn{2}-dimensional vector \eqn{\bm{b}^{ncb}} of positive numbers.
#' @param b_cb A \eqn{2}-dimensional vector \eqn{\bm{b}^{cb}} of positive numbers.
#' @param TOPICS An integer number of topics \eqn{T}.
#' @param vocab A \eqn{V}-dimensional vector of vocabulary words.
#' @param iterations An integer number of iterations. Default is 300.
#' @param seed Seed. Default is 28.
#' @param result_folder A string specifying the folder in which the results will be saved.
#'
# @seealso \code{\link{postproc_TwitterLDA}}, \code{\link{pred_TwitterLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
CGS_clickbaitLDA <- function(w, x, alphastar, betaV, beta_ncb, beta_cb,
                             b_doc, b_ncb, b_cb, TOPICS, vocab = NULL,
                             iterations = 300, seed = 28, result_folder) {
  # -------------------------------------------------------------------------- #
  # Argomenti della funzione:
  #          w : matrice D x Nmax   | n-ma parola del d-mo documento (1,...,V)
  #          x : vettore Dt x 1     | tipo del d-mo topic (1,...,Dt), 0<=Dt<=D
  #  alphastar : numero positivo    | parametro Dirichlet simmetrica
  #      betaV : numero positivo    | parametro Dirichlet simmetrica
  #   beta_ncb : numero positivo    | parametro Dirichlet simmetrica
  #    beta_cb : numero positivo    | parametro Dirichlet simmetrica
  #      b_doc : vettore 2x1        | parametro Beta
  #      b_ncb : vettore 2x1        | parametro Beta
  #       b_cb : vettore 2x1        | parametro Beta
  #     TOPICS : intero             | numero di topic
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
  cat(as.character(Sys.time()), " START\n\n", sep="")
  # Fisso il seme
  set.seed(seed)
  # Definisco alcune quantita' utili
  D <- nrow(w)
  V <- max(w)
  N <- apply(w, 1, function(x) sum(x > 0))
  if(is.null(x)) {
    x <- c()
    Dt <- 0
  } else {
    Dt <- length(x)  # documenti per cui la variabile x e' nota (ncb/cb)
  }
  if(Dt < D) x <- c(x, rep(2, D - Dt))
  # Creo matrice a partire da b_ncb e b_cb
  betaB <- c(beta_cb, beta_ncb)
  b_back <- rbind(b_cb, b_ncb)
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
  hyper <- list("w" = w, "x" = x, "alphastar" = alphastar,
                "betaV" = betaV, "betaB" = betaB,
                "b_doc" = b_doc, "b_back" = b_back,
                "iterations" = iterations, "seed" = seed,
                "T" = TOPICS, "D" = D, "Dt" = Dt, "V" = V, "N" = N, "vocab" = vocab)
  saveRDS(hyper, file.path(result_folder, "hyperparameters.RDS"))
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder")))
  # -------------------------------------------------------------------------- #
  rcpp_CGS_clickbaitLDA(hyper$w, hyper$x, hyper$alphastar, hyper$betaV,
                        hyper$betaB, hyper$b_doc, hyper$b_back,
                        hyper$iterations, hyper$T, hyper$D, hyper$Dt,
                        hyper$V, hyper$N, result_folder)
  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " END", sep="")
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
