// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_CGS_HashtagLDA
void rcpp_CGS_HashtagLDA(IntegerMatrix w, IntegerMatrix h, IntegerVector doc_users, NumericVector alphastar, NumericVector betaV, NumericVector betaH, NumericVector bH, int iterations, int TOPICS, int U, int D, int V, int H, IntegerVector N, IntegerVector L, int L_sum, int Lmax, double betaV_sum, double betaH_sum, std::string result_folder);
RcppExport SEXP _mdtm_rcpp_CGS_HashtagLDA(SEXP wSEXP, SEXP hSEXP, SEXP doc_usersSEXP, SEXP alphastarSEXP, SEXP betaVSEXP, SEXP betaHSEXP, SEXP bHSEXP, SEXP iterationsSEXP, SEXP TOPICSSEXP, SEXP USEXP, SEXP DSEXP, SEXP VSEXP, SEXP HSEXP, SEXP NSEXP, SEXP LSEXP, SEXP L_sumSEXP, SEXP LmaxSEXP, SEXP betaV_sumSEXP, SEXP betaH_sumSEXP, SEXP result_folderSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type h(hSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type doc_users(doc_usersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphastar(alphastarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaV(betaVSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaH(betaHSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bH(bHSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type TOPICS(TOPICSSEXP);
    Rcpp::traits::input_parameter< int >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type L_sum(L_sumSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< double >::type betaV_sum(betaV_sumSEXP);
    Rcpp::traits::input_parameter< double >::type betaH_sum(betaH_sumSEXP);
    Rcpp::traits::input_parameter< std::string >::type result_folder(result_folderSEXP);
    rcpp_CGS_HashtagLDA(w, h, doc_users, alphastar, betaV, betaH, bH, iterations, TOPICS, U, D, V, H, N, L, L_sum, Lmax, betaV_sum, betaH_sum, result_folder);
    return R_NilValue;
END_RCPP
}
// rcpp_CGS_LDA
void rcpp_CGS_LDA(IntegerMatrix w, NumericVector alpha, NumericVector betaV, int iterations, int TOPICS, int D, int V, IntegerVector N, int Nmax, double betaV_sum, std::string result_folder);
RcppExport SEXP _mdtm_rcpp_CGS_LDA(SEXP wSEXP, SEXP alphaSEXP, SEXP betaVSEXP, SEXP iterationsSEXP, SEXP TOPICSSEXP, SEXP DSEXP, SEXP VSEXP, SEXP NSEXP, SEXP NmaxSEXP, SEXP betaV_sumSEXP, SEXP result_folderSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaV(betaVSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type TOPICS(TOPICSSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type Nmax(NmaxSEXP);
    Rcpp::traits::input_parameter< double >::type betaV_sum(betaV_sumSEXP);
    Rcpp::traits::input_parameter< std::string >::type result_folder(result_folderSEXP);
    rcpp_CGS_LDA(w, alpha, betaV, iterations, TOPICS, D, V, N, Nmax, betaV_sum, result_folder);
    return R_NilValue;
END_RCPP
}
// rcpp_CGS_MicroblogLDA
void rcpp_CGS_MicroblogLDA(std::vector<NumericMatrix> w, IntegerVector doc_users, NumericVector alphastar, NumericVector alpha, std::vector<NumericVector> beta, NumericMatrix b, NumericVector bdelta, NumericVector bT, double alpha0, int iterations, int TOPICS, int K, int U, int D, IntegerVector V, IntegerMatrix N, NumericVector N_sum, NumericVector Nmax, NumericVector beta_sum, IntegerVector Dusers, std::string result_folder);
RcppExport SEXP _mdtm_rcpp_CGS_MicroblogLDA(SEXP wSEXP, SEXP doc_usersSEXP, SEXP alphastarSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP bdeltaSEXP, SEXP bTSEXP, SEXP alpha0SEXP, SEXP iterationsSEXP, SEXP TOPICSSEXP, SEXP KSEXP, SEXP USEXP, SEXP DSEXP, SEXP VSEXP, SEXP NSEXP, SEXP N_sumSEXP, SEXP NmaxSEXP, SEXP beta_sumSEXP, SEXP DusersSEXP, SEXP result_folderSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<NumericMatrix> >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type doc_users(doc_usersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphastar(alphastarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< std::vector<NumericVector> >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bdelta(bdeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bT(bTSEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type TOPICS(TOPICSSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type V(VSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N_sum(N_sumSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nmax(NmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_sum(beta_sumSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Dusers(DusersSEXP);
    Rcpp::traits::input_parameter< std::string >::type result_folder(result_folderSEXP);
    rcpp_CGS_MicroblogLDA(w, doc_users, alphastar, alpha, beta, b, bdelta, bT, alpha0, iterations, TOPICS, K, U, D, V, N, N_sum, Nmax, beta_sum, Dusers, result_folder);
    return R_NilValue;
END_RCPP
}
// rcpp_CGS_TwitterLDA
void rcpp_CGS_TwitterLDA(IntegerMatrix w, IntegerVector doc_users, NumericVector alphastar, NumericVector betaV, NumericVector bV, int iterations, int TOPICS, int U, int D, int V, IntegerVector N, int N_sum, int Nmax, double betaV_sum, std::string result_folder);
RcppExport SEXP _mdtm_rcpp_CGS_TwitterLDA(SEXP wSEXP, SEXP doc_usersSEXP, SEXP alphastarSEXP, SEXP betaVSEXP, SEXP bVSEXP, SEXP iterationsSEXP, SEXP TOPICSSEXP, SEXP USEXP, SEXP DSEXP, SEXP VSEXP, SEXP NSEXP, SEXP N_sumSEXP, SEXP NmaxSEXP, SEXP betaV_sumSEXP, SEXP result_folderSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type doc_users(doc_usersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphastar(alphastarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaV(betaVSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bV(bVSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type TOPICS(TOPICSSEXP);
    Rcpp::traits::input_parameter< int >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type N_sum(N_sumSEXP);
    Rcpp::traits::input_parameter< int >::type Nmax(NmaxSEXP);
    Rcpp::traits::input_parameter< double >::type betaV_sum(betaV_sumSEXP);
    Rcpp::traits::input_parameter< std::string >::type result_folder(result_folderSEXP);
    rcpp_CGS_TwitterLDA(w, doc_users, alphastar, betaV, bV, iterations, TOPICS, U, D, V, N, N_sum, Nmax, betaV_sum, result_folder);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mdtm_rcpp_CGS_HashtagLDA", (DL_FUNC) &_mdtm_rcpp_CGS_HashtagLDA, 20},
    {"_mdtm_rcpp_CGS_LDA", (DL_FUNC) &_mdtm_rcpp_CGS_LDA, 11},
    {"_mdtm_rcpp_CGS_MicroblogLDA", (DL_FUNC) &_mdtm_rcpp_CGS_MicroblogLDA, 21},
    {"_mdtm_rcpp_CGS_TwitterLDA", (DL_FUNC) &_mdtm_rcpp_CGS_TwitterLDA, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_mdtm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
