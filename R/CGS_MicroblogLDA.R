# ---------------------------------------------------------------------------- #

#' Microblog-LDA (CGS)
#'
#' @description
#' Implementation of the Collapsed Gibbs Sampler (CGS) for Microblog-LDA.
#'
#' @details
#' Implementation in R and C++.
#'
#' @param w A list of \eqn{K} \eqn{D\times N_{\text{max}}^{(k)}} matrices \eqn{\mathbf{w}^{(k)}} of integers in \eqn{\{1,\ldots,V^{(k)}\}}.
#' @param doc_users A \eqn{D}-dimensional vector of integers specifying the single author of each document. The number of authors is set to \code{U = max(doc_users)}.
#' @param alphastar A \eqn{T}-dimensional vector \eqn{\bm{\alpha}^*} of positive numbers. The length of this vector specifies the number of topics \eqn{T}.
#' @param alpha A \eqn{T}-dimensional vector \eqn{\bm{\alpha}} of positive numbers. The length of this vector specifies the number of topics \eqn{T}.
#' @param beta A list of \eqn{K} \eqn{V^{(k)}}-dimensional vectors \eqn{\bm{\beta}^{(k)}} of positive numbers.
#' @param b A list of \eqn{K} \eqn{2}-dimensional vectors \eqn{\bm{b}^{(k)}} of positive numbers.
#' @param bdelta A \eqn{2}-dimensional vector \eqn{\bm{b}^\delta} of positive numbers.
#' @param bT A \eqn{2}-dimensional vector \eqn{\bm{b}^T} of positive numbers.
#' @param alpha0 A positive number \eqn{\alpha_0}. Default is \eqn{10^{-7}}.
#' @param iterations An integer number of iterations. Default is 300.
#' @param seed Seed. Default is 28.
#' @param result_folder A string specifying the folder in which the results will be saved.
#'
#' @seealso \code{\link{postproc_MicroblogLDA}}, \code{\link{pred_MicroblogLDA}}.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
#' @export
CGS_MicroblogLDA <- function(w, doc_users, alphastar, alpha, beta, b, bdelta, bT,
                           alpha0 = 10^-7, iterations = 300, seed = 28, result_folder) {
  # -------------------------------------------------------------------------- #
  # Argomenti della funzione:
  #          w : lista di matrici D x Nmax_k  | n-ma descrittore del k-mo vocabolario nel d-mo documento (1,...,V_k)
  #  doc_users : vettore D x 1                | autore del d-mo topic (1,...,U)
  #  alphastar : vettore TOPICS x 1           | parametro Dirichlet sul simplesso dei topic
  #      alpha : vettore TOPICS x 1           | parametro Dirichlet sul simplesso dei topic (smoothing prior)
  #       beta : lista di vettori V_k x 1     | parametro Dirichlet sul simplesso delle parole del k-mo vocabolario
  #          b : lista di vettore 2x1         | parametro Beta (uno per vocabolario)
  #     bdelta : vettore 2x1                  | parametro Beta
  #         bT : vettore 2x1                  | parametro Beta
  #     alpha0 : reale positivo               | weak smoothing prior (10^-7)
  # iterations : intero                       | numero di stati della catena da campionare
  #       seed : intero                       | seme per rendere i risultati replicabili
  # -------------------------------------------------------------------------- #

  # CHECK INPUTS HERE

  # -------------------------------------------------------------------------- #
  cat(as.character(Sys.time()), " START\n\n", sep="")
  # Fisso il seme
  set.seed(seed)
  # Definisco alcune quantita' utili
  TOPICS <- length(alphastar)
  K <- length(w)
  U <- max(doc_users)
  D <- length(doc_users)
  V <- sapply(beta, length)
  N <- matrix(0, nrow = D, ncol= K)
  for (k in 1:K) {
    N[, k] <- apply(w[[k]], 1, function(x) sum(x > 0))
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
  hyper <- list("w" = w, "doc_users" = doc_users,
                "alphastar" = alphastar, "alpha" = alpha, "beta" = beta, "b" = b, "bdelta" = bdelta,
                "bT" = bT, "alpha0" = alpha0, "iterations" = iterations, "seed" = seed,
                "T" = TOPICS, "K" = K, "U" = U, "D" = D, "V" = V, "N" = N)
  saveRDS(hyper, file.path(result_folder, "hyperparameters.RDS"))
  # -------------------------------------------------------------------------- #
  rm(list = setdiff(ls(), c("hyper", "result_folder")))
  # -------------------------------------------------------------------------- #
  N_sum <- apply(hyper$N, 2, sum)
  Nmax <- apply(hyper$N, 2, max)
  beta_sum <- sapply(hyper$beta, sum)
  Dusers <- sapply(1:hyper$U, function(uu) sum(hyper$doc_users == uu))
  b <- simplify2array(hyper$b)
  rcpp_CGS_MicroblogLDA(hyper$w, hyper$doc_users-1, hyper$alphastar, hyper$alpha,
                        hyper$beta, b, hyper$bdelta, hyper$bT, hyper$alpha0,
                        hyper$iterations,hyper$T, hyper$K, hyper$U, hyper$D, hyper$V,
                        hyper$N, N_sum, Nmax, beta_sum, Dusers, result_folder)
  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " END", sep="")
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
