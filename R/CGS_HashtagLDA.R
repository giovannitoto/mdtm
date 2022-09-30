# ---------------------------------------------------------------------------- #

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
  cat("\n", as.character(Sys.time()), " START\n\n", sep="")
  # Fisso il seme
  set.seed(seed)
  # Definisco alcune quantita' utili
  TOPICS <- length(alphastar)
  U <- length(unique(doc_users))
  D <- nrow(w)
  V <- length(betaV)
  H <- length(betaH)
  N <- apply(w, 1, function(x) sum(x > 0))
  L <- apply(h, 1, function(x) sum(x > 0))
  L_sum <- sum(L)
  Lmax <- max(L)
  betaV_sum <- sum(betaV)
  betaH_sum <- sum(betaH)
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
  rcpp_CGS_HashtagLDA(w, h, doc_users-1, alphastar, betaV, betaH, bH,
                      iterations, TOPICS, U, D, V, H, N, L, L_sum, Lmax,
                      betaV_sum, betaH_sum, result_folder)
  # -------------------------------------------------------------------------- #
  cat("\n", as.character(Sys.time()), " END", sep="")
  # END FUNCTION
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
