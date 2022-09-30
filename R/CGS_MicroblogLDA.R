# ---------------------------------------------------------------------------- #

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
  cat("\n", as.character(Sys.time()), " START\n\n", sep="")
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
  N_sum <- apply(N, 2, sum)
  Nmax <- apply(N, 2, max)
  beta_sum <- sapply(beta, sum)
  Dusers <- sapply(1:U, function(uu) sum(doc_users==uu))
  b <- simplify2array(b)
  # -------------------------------------------------------------------------- #
  # Creo cartella in cui salvare gli stati della catena
  result_folder <- file.path(getwd(), result_folder)
  if(!dir.exists(result_folder)) {
    dir.create(result_folder)
    for (m in 1:iterations) {
      dir.create(file.path(result_folder, m))
    }
  } else {
    cat("\n'",result_folder,"' already exists: select another value for the 'result_folder' argument.\n", sep="")
    return(NULL)
  }
  # -------------------------------------------------------------------------- #
  rcpp_CGS_MicroblogLDA(w, doc_users-1, alphastar, alpha, beta, b,
                        bdelta, bT, alpha0, iterations, TOPICS, K, U, D, V,
                        N, N_sum, Nmax, beta_sum, Dusers, result_folder)
  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " END", sep="")
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
