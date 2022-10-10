// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "update_counts.h"
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void pred_single_HashtagLDA(IntegerMatrix w, IntegerMatrix h,
                            IntegerVector doc_users, NumericVector alphastar,
                            NumericVector betaV, NumericVector betaH, NumericVector bH,
                            int iterations, int TOPICS, int U, int D, int V, int H,
                            IntegerVector N, IntegerVector L, NumericMatrix WY1ZX,
                            NumericMatrix HY1ZX, NumericMatrix Zstar, double Yh1,
                            NumericVector HY0, std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, l, n, t, u;
  int num_count, den_count;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  int L_sum = sum(L);
  int Lmax = max(L);
  double betaV_sum = sum(betaV);
  double betaH_sum = sum(betaH);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericMatrix WY1ZX_new(V, TOPICS);
  NumericMatrix HY1ZX_new(H, TOPICS);
  NumericMatrix Zstar_new(U, TOPICS);
  double Yh1_new = 0.0;
  NumericVector HY0_new(H);
  // state of the chain
  IntegerVector zstar(D);
  IntegerMatrix yH(D, Lmax);
  // initialization of the first state
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_HashtagLDA(w, h, doc_users, alphastar, bH, TOPICS, D, N, L, zstar, yH,
                           WY1ZX_new, HY1ZX_new, Zstar_new, Yh1_new, HY0_new, true);
  // save first state of the chain
  saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(0) + "/zstar.RDS"));
  saveRDS(yH, Named("file", result_folder + "/" + std::to_string(0) + "/yH.RDS"));
  // END COUNTS AND FIRST STATE
  // ---------------------------------------------------------------------------
  // START COLLAPSED GIBBS SAMPLER
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Collapsed Gibbs Sampler (" << iterations << " iterations)\n";
  for (int m = 0; m < iterations; m++) {
    // -------------------------------------------------------------------------
    // START CGS ITERATION
    for (d = 0; d < D; d++) {
      // -----------------------------------------------------------------------
      // START DOCUMENT UPDATE
      u = doc_users(d);
      // -----------------------------------------------------------------------
      WY1ZX_new.fill(0);
      HY1ZX_new.fill(0);
      Zstar_new.fill(0);
      Yh1_new = 0.0;
      HY0_new.fill(0);
      Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) + 1;
      for (int n = 0; n < N(d); n++) {
        WY1ZX_new(w(d,n)-1, zstar(d)-1) = WY1ZX_new(w(d,n)-1, zstar(d)-1) + 1;
      }
      for (int l = 0; l < L(d); l++) {
        if(yH(d, l) == 1) {
          HY1ZX_new(h(d,l)-1, zstar(d)-1) = HY1ZX_new(h(d,l)-1, zstar(d)-1) + 1;
          Yh1 = Yh1 + 1;
        } else {
          HY0_new(h(d,l)-1) = HY0_new(h(d,l)-1) + 1;
        }
      }
      // -----------------------------------------------------------------------
      // START UPDATE zstar(d)
      // update counts
      Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) - 1;
      for (n = 0; n < N(d); n++) {
        WY1ZX_new(w(d,n)-1, zstar(d)-1) = WY1ZX_new(w(d,n)-1, zstar(d)-1) - 1;
      }
      for (l = 0; l < L(d); l++) {
        HY1ZX_new(h(d,l)-1, zstar(d)-1) = HY1ZX_new(h(d,l)-1, zstar(d)-1) - yH(d, l);
      }
      // full conditional probability
      NumericVector p_zstar = alphastar + Zstar(u, _) + Zstar_new(u, _);
      for (t = 0; t < TOPICS; t++) {
        num_count = 0;
        den_count = 0;
        p_zstar(t) = p_zstar(t) * (betaV(w(d,0)-1) + WY1ZX(w(d,0)-1, t) + WY1ZX_new(w(d,0)-1, t)) / (betaV_sum + sum(WY1ZX(_, t)) + sum(WY1ZX_new(_, t)));
        num_count++;
        den_count++;
        for (n = 1; n < N(d); n++) {
          num_count *= (w(d, n-1) == w(d, n));
          p_zstar(t) = p_zstar(t) * (betaV(w(d,n)-1) + WY1ZX(w(d,n)-1, t) + WY1ZX_new(w(d,n)-1, t) + num_count) / (betaV_sum + sum(WY1ZX(_, t)) + sum(WY1ZX_new(_, t)) + den_count);
          num_count++;
          den_count++;
        }
        num_count = 0;
        den_count = 0;
        if (yH(d, 0) == 1) {
          p_zstar(t) = p_zstar(t) * (betaH(h(d,0)-1) + HY1ZX(h(d,0)-1, t) + HY1ZX_new(h(d,0)-1, t)) / (betaH_sum + sum(HY1ZX(_, t)) + sum(HY1ZX_new(_, t)));
          num_count++;
          den_count++;
        }
        for (l = 1; l < L(d); l++) {
          num_count *= (h(d, l-1) == h(d, l));
          if (yH(d, l) == 1) {
            p_zstar(t) = p_zstar(t) * (betaH(h(d,l)-1) + HY1ZX(h(d,l)-1, t) + HY1ZX_new(h(d,l)-1, t) + num_count) / (betaH_sum + sum(HY1ZX(_, t)) + sum(HY1ZX_new(_, t)) + den_count);
            num_count++;
            den_count++;
          }
        }
      }
      // sample new value
      zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
      // update counts
      Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) + 1;
      for (n = 0; n < N(d); n++) {
        WY1ZX_new(w(d,n)-1, zstar(d)-1) = WY1ZX_new(w(d,n)-1, zstar(d)-1) + 1;
      }
      for (l = 0; l < L(d); l++) {
        HY1ZX_new(h(d,l)-1, zstar(d)-1) = HY1ZX_new(h(d,l)-1, zstar(d)-1) + yH(d, l);
      }
      // END UPDATE zstar(d)
      // -----------------------------------------------------------------------
      for (l = 0; l < L(d); l++) {
        // ---------------------------------------------------------------------
        // START UPDATE yH(d, l)
        // update counts
        if(yH(d, l) == 1) {
          HY1ZX_new(h(d,l)-1, zstar(d)-1) = HY1ZX_new(h(d,l)-1, zstar(d)-1) - 1;
          Yh1_new = Yh1_new - 1;
        } else {
          HY0_new(h(d,l)-1) = HY0_new(h(d,l)-1) - 1;
        }
        // full conditional probability
        double p_Yh0 = (bH(1) + L_sum-1 - Yh1 - Yh1_new) * (betaH(h(d,l)-1) + HY0(h(d,l)-1) + HY0_new(h(d,l)-1)) / (betaH_sum + sum(HY0) + sum(HY0_new));
        double p_Yh1 = (bH(0) + Yh1 + Yh1_new) * (betaH(h(d,l)-1) + HY1ZX(h(d,l)-1,zstar(d)-1) + HY1ZX_new(h(d,l)-1,zstar(d)-1)) / (betaH_sum + sum(HY1ZX(_,zstar(d)-1)) + sum(HY1ZX_new(_,zstar(d)-1)));
        // sample new value
        yH(d, l) = R::rbinom(1, p_Yh1 / (p_Yh0+p_Yh1));
        // update counts
        if(yH(d, l) == 1) {
          HY1ZX_new(h(d,l)-1, zstar(d)-1) = HY1ZX_new(h(d,l)-1, zstar(d)-1) + 1;
          Yh1_new = Yh1_new + 1;
        } else {
          HY0_new(h(d,l)-1) = HY0_new(h(d,l)-1) + 1;
        }
        // END UPDATE yH(d, l)
        // ---------------------------------------------------------------------
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zstar.RDS"));
    saveRDS(yH, Named("file", result_folder + "/" + std::to_string(m + 1) + "/yH.RDS"));
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
void pred_all_HashtagLDA(IntegerMatrix w, IntegerMatrix h,
                         IntegerVector doc_users, NumericVector alphastar,
                         NumericVector betaV, NumericVector betaH, NumericVector bH,
                         int iterations, int TOPICS, int U, int D, int V, int H,
                         IntegerVector N, IntegerVector L, NumericMatrix WY1ZX,
                         NumericMatrix HY1ZX, NumericMatrix Zstar, double Yh1,
                         NumericVector HY0, std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, l, n, t, u;
  int num_count, den_count;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  int L_sum = sum(L);
  int Lmax = max(L);
  double betaV_sum = sum(betaV);
  double betaH_sum = sum(betaH);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // state of the chain
  IntegerVector zstar(D);
  IntegerMatrix yH(D, Lmax);
  // initialization of the first state
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_HashtagLDA(w, h, doc_users, alphastar, bH, TOPICS, D, N, L,
                           zstar, yH, WY1ZX, HY1ZX, Zstar, Yh1, HY0, true);
  // save first state of the chain
  saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(0) + "/zstar.RDS"));
  saveRDS(yH, Named("file", result_folder + "/" + std::to_string(0) + "/yH.RDS"));
  // END COUNTS AND FIRST STATE
  // ---------------------------------------------------------------------------
  // START COLLAPSED GIBBS SAMPLER
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Collapsed Gibbs Sampler (" << iterations << " iterations)\n";
  for (int m = 0; m < iterations; m++) {
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
        WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) - 1;
      }
      for (l = 0; l < L(d); l++) {
        HY1ZX(h(d,l)-1, zstar(d)-1) = HY1ZX(h(d,l)-1, zstar(d)-1) - yH(d, l);
      }
      // full conditional probability
      NumericVector p_zstar = alphastar + Zstar(u, _);
      for (t = 0; t < TOPICS; t++) {
        num_count = 0;
        den_count = 0;
        p_zstar(t) = p_zstar(t) * (betaV(w(d,0)-1) + WY1ZX(w(d,0)-1, t)) / (betaV_sum + sum(WY1ZX(_, t)));
        num_count++;
        den_count++;
        for (n = 1; n < N(d); n++) {
          num_count *= (w(d, n-1) == w(d, n));
          p_zstar(t) = p_zstar(t) * (betaV(w(d,n)-1) + WY1ZX(w(d,n)-1, t) + num_count) / (betaV_sum + sum(WY1ZX(_, t)) + den_count);
          num_count++;
          den_count++;
        }
        num_count = 0;
        den_count = 0;
        if (yH(d, 0) == 1) {
          p_zstar(t) = p_zstar(t) * (betaH(h(d,0)-1) + HY1ZX(h(d,0)-1, t)) / (betaH_sum + sum(HY1ZX(_, t)));
          num_count++;
          den_count++;
        }
        for (l = 1; l < L(d); l++) {
          num_count *= (h(d, l-1) == h(d, l));
          if (yH(d, l) == 1) {
            p_zstar(t) = p_zstar(t) * (betaH(h(d,l)-1) + HY1ZX(h(d,l)-1, t) + num_count) / (betaH_sum + sum(HY1ZX(_, t)) + den_count);
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
        WY1ZX(w(d,n)-1, zstar(d)-1) = WY1ZX(w(d,n)-1, zstar(d)-1) + 1;
      }
      for (l = 0; l < L(d); l++) {
        HY1ZX(h(d,l)-1, zstar(d)-1) = HY1ZX(h(d,l)-1, zstar(d)-1) + yH(d, l);
      }
      // END UPDATE zstar(d)
      // -----------------------------------------------------------------------
      for (l = 0; l < L(d); l++) {
        // ---------------------------------------------------------------------
        // START UPDATE yH(d, l)
        // update counts
        if(yH(d, l) == 1) {
          HY1ZX(h(d,l)-1, zstar(d)-1) = HY1ZX(h(d,l)-1, zstar(d)-1) - 1;
          Yh1 = Yh1 - 1;
        } else {
          HY0(h(d,l)-1) = HY0(h(d,l)-1) - 1;
        }
        // full conditional probability
        double p_Yh0 = (bH(1) + L_sum-1 - Yh1) * (betaH(h(d,l)-1) + HY0(h(d,l)-1)) / (betaH_sum + sum(HY0));
        double p_Yh1 = (bH(0) + Yh1) * (betaH(h(d,l)-1) + HY1ZX(h(d,l)-1,zstar(d)-1)) / (betaH_sum + sum(HY1ZX(_,zstar(d)-1)));
        // sample new value
        yH(d, l) = R::rbinom(1, p_Yh1 / (p_Yh0+p_Yh1));
        // update counts
        if(yH(d, l) == 1) {
          HY1ZX(h(d,l)-1, zstar(d)-1) = HY1ZX(h(d,l)-1, zstar(d)-1) + 1;
          Yh1 = Yh1 + 1;
        } else {
          HY0(h(d,l)-1) = HY0(h(d,l)-1) + 1;
        }
        // END UPDATE yH(d, l)
        // ---------------------------------------------------------------------
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zstar.RDS"));
    saveRDS(yH, Named("file", result_folder + "/" + std::to_string(m + 1) + "/yH.RDS"));
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
