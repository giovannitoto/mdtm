# ---------------------------------------------------------------------------- #

dt_to_t <- function(dt, sparse = FALSE) {
  i <- j <- x <- c()
  for (r in 1:nrow(dt)) {
    if (sum(dt[r, ]) > 0) {
      x2 <- rep(1:ncol(dt), dt[r, ])
      x <- c(x, x2)
      i <- c(i, rep(r,length(x2)))
      j <- c(j, 1:length(x2))
      #print(c(r, length(i), length(j)))
    }
  }
  out <- Matrix::sparseMatrix(i, j, x = x, repr = "R",
                              dims = c(nrow(dt), max(Matrix::rowSums(dt))))
  if (sparse == FALSE) {
    out <- as.matrix(out)
    #out[out == 0] <- NA
  }
  return(out)
}

# ---------------------------------------------------------------------------- #

# https://stackoverflow.com/a/8189441

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ---------------------------------------------------------------------------- #
