// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "update_counts.h"
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void update_counts_LDA(IntegerMatrix w,
                       NumericVector alpha,
                       int TOPICS, int D, IntegerVector N,
                       IntegerMatrix&  zV,
                       NumericMatrix& WY1ZX, NumericMatrix& Z, bool update_state) {
  for (int d = 0; d < D; d++) {
    for (int n = 0; n < N(d); n++) {
      if (update_state) {
        zV(d, n) = sample(TOPICS, 1, true, alpha, true)(0);
      }
      // update counts
      WY1ZX(w(d,n)-1, zV(d,n)-1)++;
      Z(d, zV(d,n)-1)++;
    }
  }
}

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void update_counts_TwitterLDA(IntegerMatrix w, IntegerVector doc_users,
                              NumericVector alphastar, NumericVector bV,
                              int TOPICS, int D, IntegerVector N,
                              IntegerVector&  zstar, IntegerMatrix& yV,
                              NumericMatrix& WY1ZX, NumericMatrix& Zstar,
                              double& Yv1, NumericVector& WY0, bool update_state) {
  int u;
  for (int d = 0; d < D; d++) {
    u = doc_users(d);
    // sample
    if (update_state) {
      zstar(d) = sample(TOPICS, 1, true, alphastar, true)(0);
    }
    // update counts
    Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) + 1;
    for (int n = 0; n < N(d); n++) {
      // sample
      if (update_state) {
        yV(d, n) = R::rbinom(1, bV(0) / sum(bV));
      }
      // update counts
      if(yV(d, n) == 1) {
        WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) + 1;
        Yv1 = Yv1 + 1;
      } else {
        WY0(w(d,n)-1) = WY0(w(d,n)-1) + 1;
      }
    }
  }
}

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void update_counts_HashtagLDA(IntegerMatrix w, IntegerMatrix h, IntegerVector doc_users,
                              NumericVector alphastar, NumericVector bH,
                              int TOPICS, int D, IntegerVector N, IntegerVector L,
                              IntegerVector& zstar, IntegerMatrix& yH,
                              NumericMatrix& WY1ZX, NumericMatrix& HY1ZX,
                              NumericMatrix& Zstar, double& Yh1,
                              NumericVector& HY0, bool update_state) {
  int u;
  for (int d = 0; d < D; d++) {
    u = doc_users(d);
    // sample
    if (update_state) {
      zstar(d) = sample(TOPICS, 1, true, alphastar, true)(0);
    }
    // update counts
    Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) + 1;
    for (int n = 0; n < N(d); n++) {
      // update counts
      WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) + 1;
    }
    for (int l = 0; l < L(d); l++) {
      // sample
      if (update_state) {
        yH(d, l) = R::rbinom(1, bH(0) / sum(bH));
      }
      // update counts
      if(yH(d, l) == 1) {
        HY1ZX(h(d,l)-1, zstar(d)-1) = HY1ZX(h(d,l)-1, zstar(d)-1) + 1;
        Yh1 = Yh1 + 1;
      } else {
        HY0(h(d,l)-1) = HY0(h(d,l)-1) + 1;
      }
    }
  }
}

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void update_counts_MicroblogLDA(std::vector<NumericMatrix>& w,
                                IntegerVector doc_users,
                                NumericVector alphastar, NumericVector alpha,
                                std::vector<NumericVector>& beta,
                                NumericMatrix b, NumericVector bdelta,
                                NumericVector bT, double alpha0,
                                int TOPICS, int K, int D, IntegerMatrix N,
                                IntegerVector& x, IntegerVector& zstar,
                                NumericMatrix& lambda, std::vector<IntegerMatrix>& y,
                                std::vector<IntegerMatrix>& z,
                                NumericVector& X1, NumericMatrix& Zstar,
                                double& LAMBDA1, NumericMatrix& Z, NumericVector& Yv1,
                                std::vector<NumericMatrix>& WY1ZX,
                                std::vector<NumericVector>& WY0, bool update_state) {
  int u;
  for (int d = 0; d < D; d++) {
    u = doc_users(d);
    // sample values
    if (update_state) {
      x(d) = R::rbinom(1, bT(0) / sum(bT));
      zstar(d) = sample(TOPICS, 1, true, alphastar, true)(0);
    }
    // update counts
    X1(u) = X1(u) + x(d);
    Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) + 1;
    // sample and update
    for (int t = 0; t < TOPICS; t++) {
      if (update_state) {
        lambda(d, t) = R::rbinom(1, bdelta(0) / sum(bdelta));
      }
      LAMBDA1 += lambda(d, t);
    }
    NumericVector alpha_doc = alpha0 + alpha * lambda(d, _);
    for (int k = 0; k < K; k++) {
      for (int n = 0; n < N(d, k); n++) {
        // sample values
        if (update_state) {
          y[k](d, n) = R::rbinom(1, b(0,k) / sum(b(_,k)));
          z[k](d, n) = sample(TOPICS, 1, true, alpha_doc, true)(0);
        }
        // update counts
        Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) + 1;
        if (y[k](d, n) == 1) {
          if (x(d) == 1) {
            WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
          } else {
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + 1;
          }
          Yv1(k) = Yv1(k) + 1;
        } else{
          WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) + 1;
        }
      }
      // END - descriptor k
    }
    // END - document d
  }
}

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void update_counts_clickbaitLDA(IntegerMatrix w, IntegerVector& x,
                                double alphastar, double betaV, NumericVector betaB,
                                NumericVector b_doc, NumericMatrix b_back,
                                int TOPICS, int D, int Dt, int V, IntegerVector N,
                                IntegerVector&  zstar, IntegerMatrix& yV,
                                NumericMatrix& WY1ZX, NumericVector& Zstar, double& X1,
                                NumericVector& Yv1, NumericMatrix& WY0,
                                bool update_state) {
  // probability vectors
  NumericVector alphastar_vector(TOPICS, alphastar);
  for (int d = 0; d < D; d++) {
    // sample
    if (update_state) {
      zstar(d) = sample(TOPICS, 1, true, alphastar_vector, true)(0);
      if (d > Dt-1) {
        x(d) = R::rbinom(1, b_doc(0) / sum(b_doc));
      }
    }
    // update counts
    Zstar(zstar(d)-1) = Zstar(zstar(d)-1) + 1;
    X1 = X1 + x(d);
    for (int n = 0; n < N(d); n++) {
      // sample
      if (update_state) {
        yV(d, n) = R::rbinom(1, b_back(x(d),0) / sum(b_back(x(d),_)));
      }
      // update counts
      if(yV(d, n) == 1) {
        WY1ZX(zstar(d)-1, w(d,n)-1) = WY1ZX(zstar(d)-1, w(d,n)-1) + 1;
        Yv1(x(d)) = Yv1(x(d)) + 1;
      } else {
        WY0(x(d), w(d,n)-1) = WY0(x(d), w(d,n)-1) + 1;
      }
    }
  }
}

// -----------------------------------------------------------------------------
