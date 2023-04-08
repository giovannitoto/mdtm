// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "update_counts.h"
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void pred_single_LDA(IntegerMatrix w, NumericVector alpha, NumericVector betaV,
                     int iterations, int TOPICS, int D, int V, IntegerVector N,
                     NumericMatrix WY1ZX, std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, m, n;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  int Nmax = max(N);
  double betaV_sum = sum(betaV);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericMatrix WY1ZX_new(V, TOPICS);
  NumericMatrix Z(D, TOPICS);
  // state of the chain
  IntegerMatrix zV(D, Nmax);
  // initialization of the first state
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_LDA(w, alpha, TOPICS, D, N, zV, WY1ZX_new, Z, true);
  // save first state of the chain
  saveRDS(zV, Named("file", result_folder + "/" + std::to_string(0) + "/zV.RDS"));
  // END COUNTS AND FIRST STATE
  // ---------------------------------------------------------------------------
  // START COLLAPSED GIBBS SAMPLER
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Collapsed Gibbs Sampler (" << iterations << " iterations)\n";
  for (m = 0; m < iterations; m++) {
    for (d = 0; d < D; d++) {
      // -----------------------------------------------------------------------
      // START DOCUMENT UPDATE
      WY1ZX_new.fill(0);
      for (n = 0; n < N(d); n++) {
        WY1ZX_new(w(d,n)-1, zV(d,n)-1) = WY1ZX_new(w(d,n)-1, zV(d,n)-1) + 1;
      }
      for (n = 0; n < N(d); n++) {
        // ---------------------------------------------------------------------
        // START UPDATE zV(d, n)
        // update counts
        WY1ZX_new(w(d,n)-1, zV(d,n)-1) = WY1ZX_new(w(d,n)-1, zV(d,n)-1) - 1;
        Z(d, zV(d,n)-1) = Z(d, zV(d,n)-1) - 1;
        // full conditional probability
        NumericVector p_zV = (alpha + Z(d,_)) * (betaV(w(d,n)-1) + WY1ZX(w(d,n)-1,_) + WY1ZX_new(w(d,n)-1,_)) / (betaV_sum + colSums(WY1ZX) + colSums(WY1ZX_new));
        // sample new value
        zV(d, n) = sample(TOPICS, 1, true, p_zV, true)(0);
        // update counts
        WY1ZX_new(w(d,n)-1, zV(d,n)-1) = WY1ZX_new(w(d,n)-1, zV(d,n)-1) + 1;
        Z(d, zV(d,n)-1) = Z(d, zV(d,n)-1) + 1;
        // END UPDATE zV(d, n)
        // ---------------------------------------------------------------------
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(zV, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zV.RDS"));
    // -------------------------------------------------------------------------
    tt = std::time(nullptr);
    std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
    Rcout << mbstr << "  - iteration "<< m + 1 << "\n";
    // END CGS ITERATION
    // -------------------------------------------------------------------------
  }
  // END COLLAPSED GIBBS SAMPLER
  // ---------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void pred_all_LDA(IntegerMatrix w, NumericVector alpha, NumericVector betaV,
                  int iterations, int TOPICS, int D, int V, IntegerVector N,
                  NumericMatrix WY1ZX, std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, m, n;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  int Nmax = max(N);
  double betaV_sum = sum(betaV);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericMatrix Z(D, TOPICS);
  // state of the chain
  IntegerMatrix zV(D, Nmax);
  // initialization of the first state
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_LDA(w, alpha, TOPICS, D, N, zV, WY1ZX, Z, true);
  // save first state of the chain
  saveRDS(zV, Named("file", result_folder + "/" + std::to_string(0) + "/zV.RDS"));
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
      for (n = 0; n < N(d); n++) {
        // ---------------------------------------------------------------------
        // START UPDATE zV(d, n)
        // update counts
        WY1ZX(w(d,n)-1, zV(d,n)-1) = WY1ZX(w(d,n)-1, zV(d,n)-1) - 1;
        Z(d, zV(d,n)-1) = Z(d, zV(d,n)-1) - 1;
        // full conditional probability
        NumericVector p_zV = (alpha + Z(d,_)) * (betaV(w(d,n)-1) + WY1ZX(w(d,n)-1,_)) / (betaV_sum + colSums(WY1ZX));
        // sample new value
        zV(d, n) = sample(TOPICS, 1, true, p_zV, true)(0);
        // update counts
        WY1ZX(w(d,n)-1, zV(d,n)-1) = WY1ZX(w(d,n)-1, zV(d,n)-1) + 1;
        Z(d, zV(d,n)-1) = Z(d, zV(d,n)-1) + 1;
        // END UPDATE zV(d, n)
        // ---------------------------------------------------------------------
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(zV, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zV.RDS"));
    // -------------------------------------------------------------------------
    tt = std::time(nullptr);
    std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
    Rcout << mbstr << "  - iteration "<< m + 1 << "\n";
    // END CGS ITERATION
    // -------------------------------------------------------------------------
  }
  // END COLLAPSED GIBBS SAMPLER
  // ---------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------
