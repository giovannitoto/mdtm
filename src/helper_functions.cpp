// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include <ctime>
using namespace Rcpp;

// -----------------------------------------------------------------------------


//' Document-term matrix to term matrix
//'
//' @description
//' Let \eqn{V} be the dimension of the vocabulary of the collection, this
//' function convert a  \eqn{D\times V} document-term matrix into a
//' \eqn{D\times N_{\text{max}}} matrix of integers in \eqn{\{1,\ldots,V\}}.
//'
//' @param dt A \eqn{D\times V} document-term matrix.
//'
//' @return A \eqn{D\times N_{\text{max}}} matrix of integers in \eqn{\{1,\ldots,V\}}.
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix dt_to_t(const IntegerMatrix dt) {
  int c, d, n, v;
  int D = dt.nrow();
  int V = dt.ncol();
  IntegerVector N = rowSums(dt);
  int Nmax = max(N);
  IntegerMatrix out(D, Nmax);
  for (d = 0; d < D; d++) {
    n = 0;
    for (v = 0; v < V; v++) {
      if (dt(d, v) > 0) {
        for (c = 0; c < dt(d, v); c++) {
          out(d, n) = v + 1;
          n++;
        }
      }
    }
  }
  return out;
}

// -----------------------------------------------------------------------------

NumericVector get_unique_words(const NumericMatrix& w, const int& row) {
  NumericVector unique_words = w(row, _);
  unique_words = unique_words[unique_words > 0];
  unique_words = unique(unique_words);
  std::sort(unique_words.begin(), unique_words.end());
  return unique_words;
}
IntegerVector get_doc_words(const IntegerMatrix w, const int row) {
  IntegerVector doc_words = w(row, _);
  doc_words = doc_words[doc_words > 0];
  std::sort(doc_words.begin(), doc_words.end());
  return doc_words;
}

// -----------------------------------------------------------------------------
