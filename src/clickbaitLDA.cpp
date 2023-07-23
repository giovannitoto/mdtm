// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "update_counts.h"
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void rcpp_CGS_clickbaitLDA(IntegerMatrix w, IntegerVector x,
                double alphastar, double betaV, NumericVector betaB,
                NumericVector b_doc, NumericMatrix b_back,
                int iterations, int TOPICS, int D, int Dt, int V, IntegerVector N,
                std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, m, n, t;
  int num_count, den_count;
  int num1_count, den1_count, num0_count, den0_count;
  IntegerVector N_sum(2);
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  int Ntot = sum(N); // = N_sum(0) + N_sum(1)
  int Nmax = max(N);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericMatrix WY1ZX(TOPICS, V);
  NumericVector Zstar(TOPICS);
  double X1 = 0.0;
  NumericMatrix WY0(2, V);
  NumericVector Yv1(2);
  // state of the chain
  IntegerVector zstar(D);
  IntegerMatrix yV(D, Nmax);
  // initialization of the first state
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_clickbaitLDA(w, x, alphastar, betaV, betaB, b_doc, b_back,
                             TOPICS, D, Dt, V, N, zstar, yV, WY1ZX, Zstar, X1,
                             Yv1, WY0, true);
  // save first state of the chain
  saveRDS(x, Named("file", result_folder + "/" + std::to_string(0) + "/x.RDS"));
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
      // -----------------------------------------------------------------------
      // START UPDATE zstar(d)
      // update counts
      Zstar(zstar(d)-1) = Zstar(zstar(d)-1) - 1;
      for (n = 0; n < N(d); n++) {
        WY1ZX(zstar(d)-1, w(d,n)-1) = WY1ZX(zstar(d)-1, w(d,n)-1) - yV(d, n);
      }
      // full conditional probability
      NumericVector p_zstar = alphastar + Zstar;
      for (t = 0; t < TOPICS; t++) {
        num_count = 0;
        den_count = 0;
        if (yV(d, 0) == 1) {
          p_zstar(t) = p_zstar(t) * (betaV + WY1ZX(t, w(d,0)-1)) / (betaV * V + sum(WY1ZX(t, _)));
          num_count++;
          den_count++;
        }
        for (n = 1; n < N(d); n++) {
          num_count *= (w(d, n-1) == w(d, n));
          if (yV(d, n) == 1) {
            p_zstar(t) = p_zstar(t) * (betaV + WY1ZX(t, w(d,n)-1) + num_count) / (betaV * V + sum(WY1ZX(t, _)) + den_count);
            num_count++;
            den_count++;
          }
        }
      }
      // sample new value
      zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
      // update counts
      Zstar(zstar(d)-1) = Zstar(zstar(d)-1) + 1;
      for (n = 0; n < N(d); n++) {
        WY1ZX(zstar(d)-1, w(d,n)-1) = WY1ZX(zstar(d)-1, w(d,n)-1) + yV(d, n);
      }
      // END UPDATE zstar(d)
      // -----------------------------------------------------------------------
      // START UPDATE x(d)
      N_sum(0) = sum(N * (1-x));
      N_sum(1) = sum(N * x);
      if (d > Dt-1) {
        // update counts
        for (n = 0; n < N(d); n++) {
          if (yV(d, n) == 1) {
            Yv1(x(d)) = Yv1(x(d)) - 1;
          } else {
            WY0(x(d), w(d,n)-1) = WY0(x(d), w(d,n)-1) - 1;
          }
        }
        X1 = X1 - x(d);
        // full conditional probability
        double p_x1 = b_doc(0) + X1;
        double p_x0 = b_doc(1) + Ntot-1 - X1;
        num_count = 0; den_count = 0;  // entrambe nel numeratore (riusiamo le stesse variabili)
        num1_count = 0; den1_count = 0;
        num0_count = 0; den0_count = 0;
        if (yV(d, 0) == 1) {
          p_x1 *= (b_back(1,0) + Yv1(1)) / (sum(b_back(1, _)) + N_sum(1) - x(d)*N(d));
          p_x0 *= (b_back(0,0) + Yv1(0)) / (sum(b_back(0, _)) + N_sum(0) - (1-x(d))*N(d));
          num_count++;
          p_x1 *= (betaB(1) + WY0(1, w(d,0)-1)) / (betaB(1)*V + sum(WY0(1, _)));
          num1_count++; den1_count++;
        } else {
          p_x1 *= (b_back(1,1) + N_sum(1) - x(d)*N(d) - Yv1(1)) / (sum(b_back(1, _)) + N_sum(1) - x(d)*N(d));
          p_x0 *= (b_back(0,1) + N_sum(0) - (1-x(d))*N(d) - Yv1(0)) / (sum(b_back(0, _)) + N_sum(0) - (1-x(d))*N(d));
          den_count++;
          p_x0 *= (betaB(0) + WY0(0, w(d,0)-1)) / (betaB(1)*V + sum(WY0(0, _)));
          num0_count++; den0_count++;
        }
        for (n = 1; n < N(d); n++) {
          num1_count *= (w(d, n-1) == w(d, n));
          num0_count *= (w(d, n-1) == w(d, n));
          if (yV(d, n) == 1) {
            p_x1 *= (b_back(1,0) + Yv1(1) + num_count) / (sum(b_back(1, _)) + N_sum(1) - x(d)*N(d) + n);
            p_x0 *= (b_back(0,0) + Yv1(0) + num_count) / (sum(b_back(0, _)) + N_sum(0) - (1-x(d))*N(d) + n);
            num_count++;
            p_x1 *= (betaB(1) + WY0(1, w(d,0)-1) + num1_count) / (betaB(1)*V + sum(WY0(1, _)) + den1_count);
            num1_count++; den1_count++;
          } else {
            p_x1 *= (b_back(1,1) + N_sum(1) - x(d)*N(d) - Yv1(1) + den_count) / (sum(b_back(1, _)) + N_sum(1) - x(d)*N(d) + n);
            p_x0 *= (b_back(0,1) + N_sum(0) - (1-x(d))*N(d) - Yv1(0) + den_count) / (sum(b_back(0, _)) + N_sum(0) - (1-x(d))*N(d) + n);
            den_count++;
            p_x0 *= (betaB(0) + WY0(0, w(d,0)-1) + num0_count) / (betaB(1)*V + sum(WY0(0, _)) + den0_count);
            num0_count++; den0_count++;
          }
        }
        // sample new value
        x(d) = R::rbinom(1, p_x1 / (p_x0+p_x1));
        // update counts
        for (n = 0; n < N(d); n++) {
          if (yV(d, n) == 1) {
            Yv1(x(d)) = Yv1(x(d)) + 1;
          } else {
            WY0(x(d), w(d,n)-1) = WY0(x(d), w(d,n)-1) + 1;
          }
        }
        X1 = X1 + x(d);
        N_sum(0) = sum(N * (1-x));
        N_sum(1) = sum(N * x);
      }
      // END UPDATE x(d)
      // -----------------------------------------------------------------------
      for (n = 0; n < N(d); n++) {
        // ---------------------------------------------------------------------
        // START UPDATE yV(d, n)
        // update counts
        if(yV(d, n) == 1) {
          WY1ZX(zstar(d)-1, w(d,n)-1) = WY1ZX(zstar(d)-1, w(d,n)-1) - 1;
          Yv1(x(d)) = Yv1(x(d)) - 1;
        } else {
          WY0(x(d), w(d,n)-1) = WY0(x(d), w(d,n)-1) - 1;
        }
        // full conditional probability
        double p_Yv0 = (b_back(x(d), 1) + N_sum(x(d))-1 - Yv1(x(d))) * (betaB(x(d)) + WY0(x(d), w(d,n)-1)) / (betaB(x(d)) * V + sum(WY0(x(d), _)));
        double p_Yv1 = (b_back(x(d), 0) + Yv1(x(d))) * (betaV + WY1ZX(zstar(d)-1, w(d,n)-1)) / (betaV * V + sum(WY1ZX(zstar(d)-1, _)));
        // sample new value
        yV(d, n) = R::rbinom(1, p_Yv1 / (p_Yv0 + p_Yv1));
        // update counts
        if(yV(d, n) == 1) {
          WY1ZX(zstar(d)-1, w(d,n)-1) = WY1ZX(zstar(d)-1, w(d,n)-1) + 1;
          Yv1(x(d)) = Yv1(x(d)) + 1;
        } else {
          WY0(x(d), w(d,n)-1) = WY0(x(d), w(d,n)-1) + 1;
        }
        // END UPDATE yV(d, n)
        // ---------------------------------------------------------------------
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(x, Named("file", result_folder + "/" + std::to_string(m + 1) + "/x.RDS"));
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
