// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "update_counts.h"
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void rcpp_CGS_TwitterLDA(IntegerMatrix w, IntegerVector doc_users,
                NumericVector alphastar, NumericVector betaV, NumericVector bV,
                int iterations, int TOPICS, int U, int D, int V, IntegerVector N,
                std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, m, n, t, u;
  int num_count, den_count;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  int N_sum = sum(N);
  int Nmax = max(N);
  double betaV_sum = sum(betaV);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericMatrix WY1ZX(V, TOPICS);
  NumericMatrix Zstar(U, TOPICS);
  double Yv1 = 0.0;
  NumericVector WY0(V);
  // state of the chain
  IntegerVector zstar(D);
  IntegerMatrix yV(D, Nmax);
  // initialization of the first state
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_TwitterLDA(w, doc_users, alphastar, bV, TOPICS, D, N,
                           zstar, yV, WY1ZX, Zstar, Yv1, WY0, true);
  // save first state of the chain
  saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(0) + "/zstar.RDS"));
  saveRDS(yV, Named("file", result_folder + "/" + std::to_string(0) + "/yV.RDS"));
  // END COUNTS AND FIRST STATE
  // ---------------------------------------------------------------------------
  // START COLLAPSED GIBBS SAMPLER
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Collapsed Gibbs Sampler (" << iterations << " iterations)\n";
  for (m = 0; m < iterations; m++) {
    // -------------------------------------------------------------------------
    // START CGS ITERATION
    for (d = 0; d < D; d++) {
      // -----------------------------------------------------------------------
      // START DOCUMENT UPDATE
      u = doc_users(d);
      // -----------------------------------------------------------------------
      // START UPDATE zstar(d)
      // update counts
      Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) - 1;
      for (n = 0; n < N(d); n++) {
        WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) - yV(d, n);
      }
      // full conditional probability
      NumericVector p_zstar = alphastar + Zstar(u, _);
      for (t = 0; t < TOPICS; t++) {
        num_count = 0;
        den_count = 0;
        if (yV(d, 0) == 1) {
          p_zstar(t) = p_zstar(t) * (betaV(w(d,0)-1) + WY1ZX(w(d,0)-1, t)) / (betaV_sum + sum(WY1ZX(_, t)));
          num_count++;
          den_count++;
        }
        for (n = 1; n < N(d); n++) {
          num_count *= (w(d, n-1) == w(d, n));
          if (yV(d, n) == 1) {
            p_zstar(t) = p_zstar(t) * (betaV(w(d,n)-1) + WY1ZX(w(d,n)-1, t) + num_count) / (betaV_sum + sum(WY1ZX(_, t)) + den_count);
            num_count++;
            den_count++;
          }
        }
      }
      // sample new value
      zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
      // update counts
      Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) + 1;
      for (n = 0; n < N(d); n++) {
        WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) + yV(d, n);
      }
      // END UPDATE zstar(d)
      // -----------------------------------------------------------------------
      for (n = 0; n < N(d); n++) {
        // ---------------------------------------------------------------------
        // START UPDATE yV(d, n)
        // update counts
        if(yV(d, n) == 1) {
          WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) - 1;
          Yv1 = Yv1 - 1;
        } else {
          WY0(w(d,n)-1) = WY0(w(d,n)-1) - 1;
        }
        // full conditional probability
        double p_Yv0 = (bV(1) + N_sum-1 - Yv1) * (betaV(w(d,n)-1) + WY0(w(d,n)-1)) / (betaV_sum + sum(WY0));
        double p_Yv1 = (bV(0) + Yv1) * (betaV(w(d,n)-1) + WY1ZX(w(d,n)-1,zstar(d)-1)) / (betaV_sum + sum(WY1ZX(_,zstar(d)-1)));
        // sample new value
        yV(d, n) = R::rbinom(1, p_Yv1 / (p_Yv0+p_Yv1));
        // update counts
        if(yV(d, n) == 1) {
          WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) + 1;
          Yv1 = Yv1 + 1;
        } else {
          WY0(w(d,n)-1) = WY0(w(d,n)-1) + 1;
        }
        // END UPDATE yV(d, n)
        // ---------------------------------------------------------------------
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zstar.RDS"));
    saveRDS(yV, Named("file", result_folder + "/" + std::to_string(m + 1) + "/yV.RDS"));
    // -------------------------------------------------------------------------
    tt = std::time(nullptr);
    std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
    Rcout << mbstr << "  - iteration "<< m + 1 << "\n";
    // END CGS ITERATION
    // -------------------------------------------------------------------------
  }
  // // END COLLAPSED GIBBS SAMPLER
  // ---------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------
