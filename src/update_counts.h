#ifndef update_counts_H
#define update_counts_H

// -----------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

void update_counts_LDA(IntegerMatrix w,
                       NumericVector alpha,
                       int TOPICS, int D, IntegerVector N,
                       IntegerMatrix&  zV,
                       NumericMatrix& WY1ZX, NumericMatrix& Z, bool update_state);

void update_counts_TwitterLDA(IntegerMatrix w, IntegerVector doc_users,
                              NumericVector alphastar, NumericVector bV,
                              int TOPICS, int D, IntegerVector N,
                              IntegerVector&  zstar, IntegerMatrix& yV,
                              NumericMatrix& WY1ZX, NumericMatrix& Zstar,
                              double& Yv1, NumericVector& WY0, bool update_state);

void update_counts_HashtagLDA(IntegerMatrix w, IntegerMatrix h, IntegerVector doc_users,
                              NumericVector alphastar, NumericVector bH,
                              int TOPICS, int D, IntegerVector N, IntegerVector L,
                              IntegerVector& zstar, IntegerMatrix& yH,
                              NumericMatrix& WY1ZX, NumericMatrix& HY1ZX,
                              NumericMatrix& Zstar, double& Yh1,
                              NumericVector& HY0, bool update_state);

void update_counts_MicroblogLDA(std::vector<NumericMatrix> w,
                                IntegerVector doc_users, IntegerVector Dusers,
                                NumericVector alphastar, NumericVector alpha,
                                std::vector<NumericVector> beta,
                                NumericMatrix b, NumericVector bdelta,
                                NumericVector bT, double alpha0,
                                int TOPICS, int K, int U, int D, IntegerMatrix N,
                                IntegerVector& x, IntegerVector& zstar,
                                NumericMatrix& lambda, std::vector<IntegerMatrix>& y,
                                std::vector<IntegerMatrix>& z,
                                NumericVector& X1, NumericMatrix& Zstar,
                                double& LAMBDA1, NumericMatrix& Z, NumericVector& Yv1,
                                std::vector<NumericMatrix>& WY1ZX,
                                std::vector<NumericVector>& WY0, bool update_state);

// -----------------------------------------------------------------------------

#endif
