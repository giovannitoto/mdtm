// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
#include <ctime>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// -----------------------------------------------------------------------------

// NumericVector get_unique_words(const NumericMatrix& w, const int& row) {
//   NumericVector unique_words = w(row, _);
//   unique_words = unique_words[unique_words > 0];
//   unique_words = unique(unique_words);
//   std::sort(unique_words.begin(), unique_words.end());
//   return unique_words;
// }
IntegerVector get_doc_words(const IntegerMatrix w, const int row) {
  IntegerVector doc_words = w(row, _);
  doc_words = doc_words[doc_words > 0];
  std::sort(doc_words.begin(), doc_words.end());
  return doc_words;
}

// -----------------------------------------------------------------------------
