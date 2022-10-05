# ---------------------------------------------------------------------------- #
# argument order: hyperparameters, parameters, state, counts
# ---------------------------------------------------------------------------- #

loglik_LDA <- function(alpha, betaV, theta, phi, zV, WY1ZX, Z) {
  sum(t(alpha + t(Z) - 1) * log(theta)) + sum(t(betaV + WY1ZX - 1) * log(phi))
}

# ---------------------------------------------------------------------------- #

loglik_TwitterLDA <- function(alphastar, betaV, thetastar, phi, phiB, piV,
                              bV, zstar, yV, Zstar, WY1ZX, Yv1, WY0, N) {
  (bV[1] + Yv1 - 1) * log(piV) + (bV[2] + N - Yv1 -1) * log(1 - piV) +
    sum(t(alphastar + t(Zstar) - 1) * log(thetastar)) +
    sum((betaV + WY0 - 1) * log(phiB)) + sum(t(betaV + WY1ZX - 1) * log(phi))
}

# ---------------------------------------------------------------------------- #

loglik_HashtagLDA <- function(zstar, yH, thetastar, phi, psi, psiB, piH, alphastar, betaV, betaH, bH) {
  (bH[1] + Yh1 -1) * log(piH) + (bH[2] + L - Yh1 -1) * log(1-piH) +
    sum(t(alphastar + t(Zstar) - 1) * log(thetastar)) +
    sum(t(betaV + WY1ZX - 1) * log(phi)) +
    sum((betaH + HY0 - 1) * log(psiB)) + sum(t(betaH + HY1ZX - 1) * log(psi))
}

# ---------------------------------------------------------------------------- #

loglik_MicroblogLDA <- function(x, zstar, lambda, zV, yV, zH, yH,
                                thetastar, theta, phi, psi, phiB, psiB, delta, piV, piH, piT,
                                alphastar, alpha, alpha0, betaV, betaH, bdelta, bV, bH, bT) {
  (bdelta[1] + LAMBDA1 -1) * log(delta) + (bdelta[1] + D*TOPICS - LAMBDA1 -1) * log(1-delta) +
    sum((bT[1] + X1 -1) * log(piT)) + sum((bT[2] - Dusers - X1 - 1) * log(1-piT)) +
    (bV[1] + Yv1 -1) * log(piV) + (bV[2] + N - Yv1 -1) * log(1-piV) +
    (bH[1] + Yh1 -1) * log(piH) + (bH[2] + L - Yh1 -1) * log(1-piH) +
    sum(t(alphastar + t(Zstar) - 1) * log(thetastar)) +
    sum(t(alpha0 + t(lambda)*alpha + t(Z) - 1) * log(theta)) +
    sum(lgamma(TOPICS*alpha0 + apply(t(lambda) * alpha, 2, sum))) -
    sum(lgamma(alpha0 + t(lambda) * alpha)) +
    sum((betaV + WY0 - 1) * log(phiB)) + sum(t(betaV + WY1ZX - 1) * log(phi)) +
    sum((betaH + HY0 - 1) * log(psiB)) + sum(t(betaH + HY1ZX - 1) * log(psi))
}

# ---------------------------------------------------------------------------- #
