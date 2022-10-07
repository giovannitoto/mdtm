# ---------------------------------------------------------------------------- #
# argument order: hyperparameters, parameters, state, counts
# ---------------------------------------------------------------------------- #

loglik_LDA <- function(alpha, betaV, theta, phi, zV, WY1ZX, Z) {
  sum(t(alpha + t(Z) - 1) * log(theta)) + sum(t(betaV + WY1ZX - 1) * log(phi))
}

# ---------------------------------------------------------------------------- #

loglik_TwitterLDA <- function(alphastar, betaV, bV, thetastar, phi, phiB, piV,
                              zstar, yV, Zstar, WY1ZX, Yv1, WY0, N) {
  (bV[1] + Yv1 - 1) * log(piV) + (bV[2] + N - Yv1 - 1) * log(1 - piV) +
    sum(t(alphastar + t(Zstar) - 1) * log(thetastar)) +
    sum((betaV + WY0 - 1) * log(phiB)) + sum(t(betaV + WY1ZX - 1) * log(phi))
}

# ---------------------------------------------------------------------------- #

loglik_HashtagLDA <- function(alphastar, betaV, betaH, bH, thetastar, phi,
                              psi, psiB, piH, zstar, yH, WY1ZX, HY1ZX, Zstar,
                              Yh1, HY0, L) {
  (bH[1] + Yh1 - 1) * log(piH) + (bH[2] + L - Yh1 - 1) * log(1 - piH) +
    sum(t(alphastar + t(Zstar) - 1) * log(thetastar)) +
    sum(t(betaV + WY1ZX - 1) * log(phi)) +
    sum((betaH + HY0 - 1) * log(psiB)) + sum(t(betaH + HY1ZX - 1) * log(psi))
}

# ---------------------------------------------------------------------------- #

loglik_MicroblogLDA <- function(alphastar, alpha, alpha0, beta, bdelta, b, bT,
                                thetastar, theta, phi, phiB, delta, Pi, piT,
                                x, zstar, lambda, y, z,
                                X1, Zstar, LAMBDA1, Z, Yv1, WY1ZX, WY0,
                                D, TOPICS, Dusers, N, K) {
  out <- (bdelta[1] + LAMBDA1 - 1) * log(delta) + (bdelta[1] + D*TOPICS - LAMBDA1 - 1) * log(1 - delta) +
    sum((bT[1] + X1 - 1) * log(piT)) + sum((bT[2] + Dusers - X1 - 1) * log(1 - piT)) +
    sum(t(alphastar + t(Zstar) - 1) * log(thetastar)) +
    sum(t(alpha0 + t(lambda)*alpha + t(Z) - 1) * log(theta)) +
    sum(lgamma(TOPICS*alpha0 + apply(t(lambda) * alpha, 2, sum))) -
    sum(lgamma(alpha0 + t(lambda) * alpha))
    for (k in 1:K) {
      out <- out + sum((beta[[k]] + WY0[[k]] - 1) * log(phiB[[k]])) +
        sum(t(beta[[k]] + WY1ZX[[k]] - 1) * log(phi[[k]])) +
        (b[[k]][1] + Yv1[k] -1) * log(Pi[[k]]) + (b[[k]][2] + N[k] - Yv1[k] - 1) * log(1 - Pi[[k]])
    }
    return(out)
}

# ---------------------------------------------------------------------------- #
