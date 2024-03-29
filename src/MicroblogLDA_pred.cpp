// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "update_counts.h"
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void pred_single_MicroblogLDA(std::vector<NumericMatrix>& w,
                              IntegerVector doc_users,
                              NumericVector alphastar, NumericVector alpha,
                              std::vector<NumericVector>& beta,
                              NumericMatrix b, NumericVector bdelta,
                              NumericVector bT, double alpha0, int iterations,
                              int TOPICS, int K, int U, int D, IntegerVector V,
                              IntegerMatrix N, NumericVector N_sum,
                              IntegerVector Dusers,
                              NumericVector& X1, NumericMatrix& Zstar,
                              double& LAMBDA1, NumericVector& Yv1,
                              std::vector<NumericMatrix>& WY1ZX,
                              std::vector<NumericVector>& WY0,
                              std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, k, n, t, u;
  int num_count, den_count; //, num0_count, den0_count;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  NumericVector Nmax(K);
  NumericVector beta_sum(K);
  for (k = 0; k < K; k++) {
    Nmax(k) = max(N(_, k));
    beta_sum(k) = sum(beta[k]);
  }
  int D_tot = sum(Dusers);
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericVector X1_new(U);
  NumericMatrix Zstar_new(U, TOPICS);
  double LAMBDA1_new = 0;
  NumericMatrix Z_new(D, TOPICS);
  NumericVector Yv1_new(K);
  std::vector<NumericMatrix> WY1ZX_new(K);
  std::vector<NumericVector> WY0_new(K);
  for (k = 0; k < K; k++) {
    WY1ZX_new[k] = NumericMatrix(V(k), TOPICS);
    WY0_new[k] = NumericVector(V(k));
  }
  // state of the chain
  IntegerVector x(D);
  IntegerVector zstar(D);
  NumericMatrix lambda(D, TOPICS);
  std::vector<IntegerMatrix> y(K);
  std::vector<IntegerMatrix> z(K);
  for (k = 0; k < K; k++) {
    y[k] = IntegerMatrix(D, Nmax(k));
    z[k] = IntegerMatrix(D, Nmax(k));
  }
  // sample first state of the chain
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_MicroblogLDA(w, doc_users, alphastar,
                             alpha, beta, b, bdelta, bT, alpha0,
                             TOPICS, K, D, N,
                             x, zstar, lambda, y, z,
                             X1_new, Zstar_new, LAMBDA1_new, Z_new,
                             Yv1_new, WY1ZX_new, WY0_new, true);
  // save first state of the chain
  saveRDS(x, Named("file", result_folder + "/" + std::to_string(0) + "/x.RDS"));
  saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(0) + "/zstar.RDS"));
  saveRDS(lambda, Named("file", result_folder + "/" + std::to_string(0) + "/lambda.RDS"));
  for (k = 0; k < K; k++) {
    saveRDS(y[k], Named("file", result_folder + "/" + std::to_string(0) + "/y" + std::to_string(k + 1) + ".RDS"));
    saveRDS(z[k], Named("file", result_folder + "/" + std::to_string(0) + "/z" + std::to_string(k + 1) + ".RDS"));
  }
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
      X1_new.fill(0);
      Zstar_new.fill(0);
      LAMBDA1_new = 0;
      Yv1_new.fill(0);
      for (k = 0; k < K; k++) {
        WY1ZX_new[k].fill(0);
        WY0_new[k].fill(0);
      }
      X1_new(u) = X1_new(u) + x(d);
      Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) + 1;
      for (t = 0; t < TOPICS; t++) {
        LAMBDA1_new += lambda(d, t);
      }
      for (k = 0; k < K; k++) {
        for (n = 0; n < N(d, k); n++) {
          if (y[k](d, n) == 1) {
            if (x(d) == 1) {
              WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
            } else {
              WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) + 1;
            }
            Yv1_new(k) = Yv1_new(k) + 1;
          } else{
            WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) + 1;
          }
        }
      }
      // -----------------------------------------------------------------------
      // START UPDATE x(d)
      // update counts
      if (x(d) == 1) {
          X1_new(u) = X1_new(u) - 1;
          for (k = 0; k < K; k++) {
            for (n = 0; n < N(d, k); n++) {
              WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) - y[k](d, n);
            }
          }
      } else {
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) - y[k](d, n);
          }
        }
      }
      // full conditional probability
      double p_x1 = bT(0) + X1(u) + X1_new(u);
      double p_x0 = bT(1) + Dusers(u)-1 - X1(u) - X1_new(u);
      for (k = 0; k < K; k++) {
        NumericVector num1_count(TOPICS);
        NumericVector den1_count(TOPICS);
        for (t = 0; t < TOPICS; t++) {
          if (y[k](d, 0) == 1) {
            if (z[k](d, 0)-1 == t) {
              p_x1 *= (beta[k](w[k](d,0)-1) + WY1ZX[k](w[k](d,0)-1,t) + WY1ZX_new[k](w[k](d,0)-1,t)) / (beta_sum(k) + sum(WY1ZX[k](_, t)) + sum(WY1ZX_new[k](_, t)));
              num1_count(t) = num1_count(t) + 1;
              den1_count(t) = den1_count(t) + 1;
            }
          }
          for (n = 1; n < N(d, k); n++) {
            if (w[k](d, n-1) != w[k](d, n)) {
              num1_count.fill(0);
            }
            if (y[k](d, n) == 1) {
              if (z[k](d, n)-1 == t) {
                p_x1 *= (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1,t) + WY1ZX_new[k](w[k](d,n)-1,t) + num1_count(t)) / (beta_sum(k) + sum(WY1ZX[k](_, t)) + sum(WY1ZX_new[k](_, t)) + den1_count(t));
                num1_count(t) = num1_count(t) + 1;
                den1_count(t) = den1_count(t) + 1;
              }
            }
          }
        }
        int num0_count = 0;
        int den0_count = 0;
        if (y[k](d, 0) == 1) {
          p_x0 *= (beta[k](w[k](d,0)-1) + WY1ZX[k](w[k](d,0)-1,zstar(d)-1) + WY1ZX_new[k](w[k](d,0)-1,zstar(d)-1)) / (beta_sum(k) + sum(WY1ZX[k](_,zstar(d)-1)) + sum(WY1ZX_new[k](_,zstar(d)-1)));
          num0_count++;
          den0_count++;
        }
        for (n = 1; n < N(d, k); n++) {
          num0_count *= (w[k](d, n-1) == w[k](d, n));
          if (y[k](d, n) == 1) {
            p_x0 *= (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1,zstar(d)-1) + WY1ZX_new[k](w[k](d,n)-1,zstar(d)-1) + num0_count) / (beta_sum(k) + sum(WY1ZX[k](_,zstar(d)-1)) + sum(WY1ZX_new[k](_,zstar(d)-1)) + den0_count);
            num0_count++;
            den0_count++;
          }
        }
        // Rcout << p_x0 << "\t" << p_x1 << "  (" << k+1 << ")\n";
      }
      // sample new value
      x(d) = R::rbinom(1, p_x1 / (p_x0+p_x1));
      // update counts
      if (x(d) == 1) {
        X1_new(u) = X1_new(u) + 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) + y[k](d, n);
          }
        }
      } else {
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) + y[k](d, n);
          }
        }
      }
      // END UPDATE x(d)
      // -----------------------------------------------------------------------
      // START UPDATE zstar(d)
      if (x(d) == 1) {
        // update counts
        Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) - 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) - y[k](d, n);
          }
        }
        // full conditional probability
        NumericVector p_zstar = alphastar + Zstar(u, _) + Zstar_new(u, _);
        // sample new value
        zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
        // update counts
        Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) + 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) + y[k](d, n);
          }
        }
      } else {
        // update counts
        Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) - 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) - y[k](d, n);
          }
        }
        // full conditional probability
        NumericVector p_zstar = alphastar + Zstar(u, _) + Zstar_new(u, _);
        for (k = 0; k < K; k++) {
          for (t = 0; t < TOPICS; t++) {
            num_count = 0;
            den_count = 0;
            if (y[k](d, 0) == 1) {
              p_zstar(t) = p_zstar(t) * (beta[k](w[k](d,0)-1) + WY1ZX[k](w[k](d,0)-1, t) + WY1ZX_new[k](w[k](d,0)-1, t)) / (beta_sum(k) + sum(WY1ZX[k](_, t)) + sum(WY1ZX_new[k](_, t)));
              num_count++;
              den_count++;
            }
            for (n = 1; n < N(d, k); n++) {
              num_count *= (w[k](d, n-1) == w[k](d, n));
              if (y[k](d, n) == 1) {
                p_zstar(t) = p_zstar(t) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, t) + WY1ZX_new[k](w[k](d,n)-1, t) + num_count) / (beta_sum(k) + sum(WY1ZX[k](_, t)) + sum(WY1ZX_new[k](_, t)) + den_count);
                num_count++;
                den_count++;
              }
            }
          }
        }
        // sample new value
        zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
        // update counts
        Zstar_new(u, zstar(d)-1) = Zstar_new(u, zstar(d)-1) + 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) + y[k](d, n);
          }
        }
      }
      // END UPDATE zstar(d)
      // -----------------------------------------------------------------------
      // START UPDATE lambda(d, _)
      for (t = 0; t < TOPICS; t++) {
        // update counts
        LAMBDA1_new -= lambda(d, t);
        // full conditional probability
        double p_lambda1 = bdelta(0) + LAMBDA1 + LAMBDA1_new;
        double p_lambda0 = bdelta(1) + D_tot*TOPICS-1 - LAMBDA1 - LAMBDA1_new;
        double lambda_alpha1 = sum(alpha * lambda(d, _));
        double lambda_alpha0 = lambda_alpha1 - alpha(t) * lambda(d, t);
        for (num_count = 0; num_count < Z_new(d,t); num_count++) {
          p_lambda1 *= (alpha0 + alpha(t) + num_count) / (TOPICS*alpha0 + lambda_alpha1 + num_count);
          p_lambda0 *= (alpha0 + num_count) / (TOPICS*alpha0 + lambda_alpha0 + num_count);
        }
        for (num_count = Z_new(d,t); num_count < sum(N(d, _)); num_count++) {
          p_lambda1 /= (TOPICS*alpha0 + lambda_alpha1 + num_count);
          p_lambda0 /= (TOPICS*alpha0 + lambda_alpha0 + num_count);
        }
        // sample new value
        lambda(d, t) = R::rbinom(1, p_lambda1 / (p_lambda0+p_lambda1));
        // update counts
        LAMBDA1_new += lambda(d, t);
      }
      // END UPDATE lambda(d, _)
      // -----------------------------------------------------------------------
      for (k = 0; k < K; k++) {
        for (n = 0; n < N(d, k); n++) {
          // ---------------------------------------------------------------------
          // START UPDATE y[k](d, n)
          if (y[k](d, n) == 1) {
            if (x(d) == 1) {
              // update counts
              WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) - 1;
              Yv1_new(k) = Yv1_new(k) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k) - Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1) + WY0_new[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]) + sum(WY0_new[k]));
              double p_y1 = (b(0,k) + Yv1(k) + Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, z[k](d,n)-1)) + sum(WY1ZX_new[k](_, z[k](d,n)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
                Yv1_new(k) = Yv1_new(k) + 1;
              } else {
                WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) + 1;
              }
            } else {
              // update counts
              WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) - 1;
              Yv1_new(k) = Yv1_new(k) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k) - Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1) + WY0_new[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]) + sum(WY0_new[k]));
              double p_y1 = (b(0,k) + Yv1(k) + Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, zstar(d)-1)) + sum(WY1ZX_new[k](_, zstar(d)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) + 1;
                Yv1_new(k) = Yv1_new(k) + 1;
              } else {
                WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) + 1;
              }
            }
          } else {
            if (x(d) == 1) {
              // update counts
              WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k) - Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1) + WY0_new[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]) + sum(WY0_new[k]));
              double p_y1 = (b(0,k) + Yv1(k) + Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, z[k](d,n)-1)) + sum(WY1ZX_new[k](_, z[k](d,n)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
                Yv1_new(k) = Yv1_new(k) + 1;
              } else {
                WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) + 1;
              }
            } else {
              // update counts
              WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k) - Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1) + WY0_new[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]) + sum(WY0_new[k]));
              double p_y1 = (b(0,k) + Yv1(k) + Yv1_new(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, zstar(d)-1)) + sum(WY1ZX_new[k](_, zstar(d)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX_new[k](w[k](d,n)-1, zstar(d)-1) + 1;
                Yv1_new(k) = Yv1_new(k) + 1;
              } else {
                WY0_new[k](w[k](d,n)-1) = WY0_new[k](w[k](d,n)-1) + 1;
              }
            }
          }
          // END UPDATE y[k](d, n)
          // ---------------------------------------------------------------------
          // START UPDATE z[k](d, n)
          if (y[k](d, n) == 0) {
            // update counts
            Z_new(d, z[k](d,n)-1) = Z_new(d, z[k](d,n)-1) - 1;
            // full conditional probability
            NumericVector p_z = alpha0 + lambda(d, _) * alpha + Z_new(d, _);
            // sample new value
            z[k](d, n) = sample(TOPICS, 1, true, p_z, true)(0);
            // update counts
            Z_new(d, z[k](d,n)-1) = Z_new(d, z[k](d,n)-1) + 1;
          } else if (x(d) == 1) {
            // update counts
            WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) - 1;
            Z_new(d, z[k](d,n)-1) = Z_new(d, z[k](d,n)-1) - 1;
            // full conditional probability
            NumericVector p_z = (alpha0 + lambda(d, _) * alpha + Z_new(d, _)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1,_) + WY1ZX_new[k](w[k](d,n)-1,_)) / (beta_sum(k) + colSums(WY1ZX[k]) + colSums(WY1ZX_new[k]));
            // sample new value
            z[k](d, n) = sample(TOPICS, 1, true, p_z, true)(0);
            // update counts
            WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX_new[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
            Z_new(d, z[k](d,n)-1) = Z_new(d, z[k](d,n)-1) + 1;
          } else {
            // update counts
            Z_new(d, z[k](d,n)-1) = Z_new(d, z[k](d,n)-1) - 1;
            // full conditional probability
            NumericVector p_z = alpha0 + lambda(d, _) * alpha + Z_new(d, _);
            // sample new value
            z[k](d, n) = sample(TOPICS, 1, true, p_z, true)(0);
            // update counts
            Z_new(d, z[k](d,n)-1) = Z_new(d, z[k](d,n)-1) + 1;
          }
          // END UPDATE z[k](d, n)
          // ---------------------------------------------------------------------
        }
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(x, Named("file", result_folder + "/" + std::to_string(m + 1) + "/x.RDS"));
    saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zstar.RDS"));
    saveRDS(lambda, Named("file", result_folder + "/" + std::to_string(m + 1) + "/lambda.RDS"));
    for (k = 0; k < K; k++) {
      saveRDS(y[k], Named("file", result_folder + "/" + std::to_string(m + 1) + "/y" + std::to_string(k + 1) + ".RDS"));
      saveRDS(z[k], Named("file", result_folder + "/" + std::to_string(m + 1) + "/z" + std::to_string(k + 1) + ".RDS"));
    }
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
void pred_all_MicroblogLDA(std::vector<NumericMatrix>& w,
                           IntegerVector doc_users,
                           NumericVector alphastar, NumericVector alpha,
                           std::vector<NumericVector>& beta,
                           NumericMatrix b, NumericVector bdelta,
                           NumericVector bT, double alpha0, int iterations,
                           int TOPICS, int K, int U, int D, IntegerVector V,
                           IntegerMatrix N, NumericVector N_sum,
                           IntegerVector Dusers,
                           NumericVector& X1, NumericMatrix& Zstar,
                           double& LAMBDA1, NumericVector& Yv1,
                           std::vector<NumericMatrix>& WY1ZX,
                           std::vector<NumericVector>& WY0,
                           std::string result_folder) {
  // ---------------------------------------------------------------------------
  // import saveRDS
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];
  // init
  int d, k, n, t, u;
  int num_count, den_count; //, num0_count, den0_count;
  std::time_t tt;
  char mbstr[100];
  // ---------------------------------------------------------------------------
  NumericVector Nmax(K);
  NumericVector beta_sum(K);
  for (k = 0; k < K; k++) {
    Nmax(k) = max(N(_, k));
    beta_sum(k) = sum(beta[k]);
  }
  // ---------------------------------------------------------------------------
  // START COUNTS AND FIRST STATE
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Creation of the matrices of counts\n";
  // matrices of counts
  NumericMatrix Z(D, TOPICS);
  // state of the chain
  IntegerVector x(D);
  IntegerVector zstar(D);
  NumericMatrix lambda(D, TOPICS);
  std::vector<IntegerMatrix> y(K);
  std::vector<IntegerMatrix> z(K);
  for (k = 0; k < K; k++) {
    y[k] = IntegerMatrix(D, Nmax(k));
    z[k] = IntegerMatrix(D, Nmax(k));
  }
  // sample first state of the chain
  tt = std::time(nullptr);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M:%S", std::localtime(&tt));
  Rcout << mbstr << " Generation of the first state\n";
  update_counts_MicroblogLDA(w, doc_users, alphastar,
                             alpha, beta, b, bdelta, bT, alpha0,
                             TOPICS, K, D, N,
                             x, zstar, lambda, y, z,
                             X1, Zstar, LAMBDA1, Z, Yv1, WY1ZX, WY0, true);
  // save first state of the chain
  saveRDS(x, Named("file", result_folder + "/" + std::to_string(0) + "/x.RDS"));
  saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(0) + "/zstar.RDS"));
  saveRDS(lambda, Named("file", result_folder + "/" + std::to_string(0) + "/lambda.RDS"));
  for (k = 0; k < K; k++) {
    saveRDS(y[k], Named("file", result_folder + "/" + std::to_string(0) + "/y" + std::to_string(k + 1) + ".RDS"));
    saveRDS(z[k], Named("file", result_folder + "/" + std::to_string(0) + "/z" + std::to_string(k + 1) + ".RDS"));
  }
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
      // START UPDATE x(d)
      // update counts
      if (x(d) == 1) {
          X1(u) = X1(u) - 1;
          for (k = 0; k < K; k++) {
            for (n = 0; n < N(d, k); n++) {
              WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) - y[k](d, n);
            }
          }
      } else {
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) - y[k](d, n);
          }
        }
      }
      // full conditional probability
      double p_x1 = bT(0) + X1(u);
      double p_x0 = bT(1) + Dusers(u)-1 - X1(u);
      for (k = 0; k < K; k++) {
        NumericVector num1_count(TOPICS);
        NumericVector den1_count(TOPICS);
        for (t = 0; t < TOPICS; t++) {
          if (y[k](d, 0) == 1) {
            if (z[k](d, 0)-1 == t) {
              p_x1 *= (beta[k](w[k](d,0)-1) + WY1ZX[k](w[k](d,0)-1,t)) / (beta_sum(k) + sum(WY1ZX[k](_, t)));
              num1_count(t) = num1_count(t) + 1;
              den1_count(t) = den1_count(t) + 1;
            }
          }
          for (n = 1; n < N(d, k); n++) {
            if (w[k](d, n-1) != w[k](d, n)) {
              num1_count.fill(0);
            }
            if (y[k](d, n) == 1) {
              if (z[k](d, n)-1 == t) {
                p_x1 *= (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1,t) + num1_count(t)) / (beta_sum(k) + sum(WY1ZX[k](_, t)) + den1_count(t));
                num1_count(t) = num1_count(t) + 1;
                den1_count(t) = den1_count(t) + 1;
              }
            }
          }
        }
        int num0_count = 0;
        int den0_count = 0;
        if (y[k](d, 0) == 1) {
          p_x0 *= (beta[k](w[k](d,0)-1) + WY1ZX[k](w[k](d,0)-1,zstar(d)-1)) / (beta_sum(k) + sum(WY1ZX[k](_,zstar(d)-1)));
          num0_count++;
          den0_count++;
        }
        for (n = 1; n < N(d, k); n++) {
          num0_count *= (w[k](d, n-1) == w[k](d, n));
          if (y[k](d, n) == 1) {
            p_x0 *= (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1,zstar(d)-1) + num0_count) / (beta_sum(k) + sum(WY1ZX[k](_,zstar(d)-1)) + den0_count);
            num0_count++;
            den0_count++;
          }
        }
        // Rcout << p_x0 << "\t" << p_x1 << "  (" << k+1 << ")\n";
      }
      // sample new value
      x(d) = R::rbinom(1, p_x1 / (p_x0+p_x1));
      // update counts
      if (x(d) == 1) {
        X1(u) = X1(u) + 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + y[k](d, n);
          }
        }
      } else {
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + y[k](d, n);
          }
        }
      }
      // END UPDATE x(d)
      // -----------------------------------------------------------------------
      // START UPDATE zstar(d)
      if (x(d) == 1) {
        // update counts
        Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) - 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) - y[k](d, n);
          }
        }
        // full conditional probability
        NumericVector p_zstar = alphastar + Zstar(u, _);
        // sample new value
        zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
        // update counts
        Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) + 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + y[k](d, n);
          }
        }
      } else {
        // update counts
        Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) - 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) - y[k](d, n);
          }
        }
        // full conditional probability
        NumericVector p_zstar = alphastar + Zstar(u, _);
        for (k = 0; k < K; k++) {
          for (t = 0; t < TOPICS; t++) {
            num_count = 0;
            den_count = 0;
            if (y[k](d, 0) == 1) {
              p_zstar(t) = p_zstar(t) * (beta[k](w[k](d,0)-1) + WY1ZX[k](w[k](d,0)-1, t)) / (beta_sum(k) + sum(WY1ZX[k](_, t)));
              num_count++;
              den_count++;
            }
            for (n = 1; n < N(d, k); n++) {
              num_count *= (w[k](d, n-1) == w[k](d, n));
              if (y[k](d, n) == 1) {
                p_zstar(t) = p_zstar(t) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, t) + num_count) / (beta_sum(k) + sum(WY1ZX[k](_, t)) + den_count);
                num_count++;
                den_count++;
              }
            }
          }
        }
        // sample new value
        zstar(d) = sample(TOPICS, 1, true, p_zstar, true)(0);
        // update counts
        Zstar(u, zstar(d)-1) = Zstar(u, zstar(d)-1) + 1;
        for (k = 0; k < K; k++) {
          for (n = 0; n < N(d, k); n++) {
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + y[k](d, n);
          }
        }
      }
      // END UPDATE zstar(d)
      // -----------------------------------------------------------------------
      // START UPDATE lambda(d, _)
      for (t = 0; t < TOPICS; t++) {
        // update counts
        LAMBDA1 -= lambda(d, t);
        // full conditional probability
        double p_lambda1 = bdelta(0) + LAMBDA1;
        double p_lambda0 = bdelta(1) + D*TOPICS-1 - LAMBDA1;
        double lambda_alpha1 = sum(alpha * lambda(d, _));
        double lambda_alpha0 = lambda_alpha1 - alpha(t) * lambda(d, t);
        for (num_count = 0; num_count < Z(d,t); num_count++) {
          p_lambda1 *= (alpha0 + alpha(t) + num_count) / (TOPICS*alpha0 + lambda_alpha1 + num_count);
          p_lambda0 *= (alpha0 + num_count) / (TOPICS*alpha0 + lambda_alpha0 + num_count);
        }
        for (num_count = Z(d,t); num_count < sum(N(d, _)); num_count++) {
          p_lambda1 /= (TOPICS*alpha0 + lambda_alpha1 + num_count);
          p_lambda0 /= (TOPICS*alpha0 + lambda_alpha0 + num_count);
        }
        // sample new value
        lambda(d, t) = R::rbinom(1, p_lambda1 / (p_lambda0+p_lambda1));
        // update counts
        LAMBDA1 += lambda(d, t);
      }
      // END UPDATE lambda(d, _)
      // -----------------------------------------------------------------------
      for (k = 0; k < K; k++) {
        for (n = 0; n < N(d, k); n++) {
          // ---------------------------------------------------------------------
          // START UPDATE y[k](d, n)
          if (y[k](d, n) == 1) {
            if (x(d) == 1) {
              // update counts
              WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) - 1;
              Yv1(k) = Yv1(k) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]));
              double p_y1 = (b(0,k) + Yv1(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, z[k](d,n)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
                Yv1(k) = Yv1(k) + 1;
              } else {
                WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) + 1;
              }
            } else {
              // update counts
              WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) - 1;
              Yv1(k) = Yv1(k) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]));
              double p_y1 = (b(0,k) + Yv1(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, zstar(d)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, zstar(d)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + 1;
                Yv1(k) = Yv1(k) + 1;
              } else {
                WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) + 1;
              }
            }
          } else {
            if (x(d) == 1) {
              // update counts
              WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]));
              double p_y1 = (b(0,k) + Yv1(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, z[k](d,n)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
                Yv1(k) = Yv1(k) + 1;
              } else {
                WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) + 1;
              }
            } else {
              // update counts
              WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) - 1;
              // full conditional probability
              double p_y0 = (b(1,k) + N_sum(k)-1 - Yv1(k)) * (beta[k](w[k](d,n)-1) + WY0[k](w[k](d,n)-1)) / (beta_sum(k) + sum(WY0[k]));
              double p_y1 = (b(0,k) + Yv1(k)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1, zstar(d)-1)) / (beta_sum(k) + sum(WY1ZX[k](_, zstar(d)-1)));
              // sample new value
              y[k](d, n) = R::rbinom(1, p_y1 / (p_y0+p_y1));
              // update counts
              if (y[k](d, n) == 1) {
                WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + 1;
                Yv1(k) = Yv1(k) + 1;
              } else {
                WY0[k](w[k](d,n)-1) = WY0[k](w[k](d,n)-1) + 1;
              }
            }
          }
          // END UPDATE y[k](d, n)
          // ---------------------------------------------------------------------
          // START UPDATE z[k](d, n)
          if (y[k](d, n) == 0) {
            // update counts
            Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) - 1;
            // full conditional probability
            NumericVector p_z = alpha0 + lambda(d, _) * alpha + Z(d, _);
            // sample new value
            z[k](d, n) = sample(TOPICS, 1, true, p_z, true)(0);
            // update counts
            Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) + 1;
          } else if (x(d) == 1) {
            // update counts
            WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) - 1;
            Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) - 1;
            // full conditional probability
            NumericVector p_z = (alpha0 + lambda(d, _) * alpha + Z(d, _)) * (beta[k](w[k](d,n)-1) + WY1ZX[k](w[k](d,n)-1,_)) / (beta_sum(k) + colSums(WY1ZX[k]));
            // sample new value
            z[k](d, n) = sample(TOPICS, 1, true, p_z, true)(0);
            // update counts
            WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) = WY1ZX[k](w[k](d,n)-1, z[k](d,n)-1) + 1;
            Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) + 1;
          } else {
            // update counts
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) - 1;
            Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) - 1;
            // full conditional probability
            NumericVector p_z = alpha0 + lambda(d, _) * alpha + Z(d, _);
            // sample new value
            z[k](d, n) = sample(TOPICS, 1, true, p_z, true)(0);
            // update counts
            WY1ZX[k](w[k](d,n)-1, zstar(d)-1) = WY1ZX[k](w[k](d,n)-1, zstar(d)-1) + 1;
            Z(d, z[k](d,n)-1) = Z(d, z[k](d,n)-1) + 1;
          }
          // END UPDATE z[k](d, n)
          // ---------------------------------------------------------------------
        }
      }
      // END DOCUMENT UPDATE
      // -----------------------------------------------------------------------
    }
    // save m-th state of the chain
    saveRDS(x, Named("file", result_folder + "/" + std::to_string(m + 1) + "/x.RDS"));
    saveRDS(zstar, Named("file", result_folder + "/" + std::to_string(m + 1) + "/zstar.RDS"));
    saveRDS(lambda, Named("file", result_folder + "/" + std::to_string(m + 1) + "/lambda.RDS"));
    for (k = 0; k < K; k++) {
      saveRDS(y[k], Named("file", result_folder + "/" + std::to_string(m + 1) + "/y" + std::to_string(k + 1) + ".RDS"));
      saveRDS(z[k], Named("file", result_folder + "/" + std::to_string(m + 1) + "/z" + std::to_string(k + 1) + ".RDS"));
    }
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
