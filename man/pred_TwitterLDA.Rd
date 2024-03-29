% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_TLDA.R
\name{pred_TwitterLDA}
\alias{pred_TwitterLDA}
\title{Prediction of new documents using Twitter-LDA (CGS)}
\usage{
pred_TwitterLDA(
  w,
  doc_users,
  betaV_new = NULL,
  postproc_file,
  single_doc = TRUE,
  iterations = 300,
  seed = 28,
  result_folder
)
}
\arguments{
\item{w}{A \eqn{D\times N_{\text{max}}} matrix \eqn{\mathbf{w}} of integers in \eqn{\{1,\ldots,V+V_{\text{new}}\}}.}

\item{doc_users}{A \eqn{D}-dimensional vector of integers specifying the single author of each document. The number of authors is set to \code{U = max(doc_users)}.}

\item{betaV_new}{A \eqn{V_{\text{new}}}-dimensional vector \eqn{\bm{\beta}_{\text{new}}^V} of positive numbers.}

\item{postproc_file}{A string specifying the RDS file in which the results of post-processing obtained using the function \code{\link{postproc_TwitterLDA}} are saved.}

\item{single_doc}{Logical: if \code{TRUE}, the algorithm predicts one document at a time, otherwise all at once. Default is \code{TRUE}.}

\item{iterations}{An integer number of iterations. Default is 300.}

\item{seed}{Seed. Default is 28.}

\item{result_folder}{A string specifying the folder in which the results will be saved.}
}
\description{
Predict new documents using an already-estimated Twitter-LDA.
}
\details{
Implementation in R and C++.
}
\note{
This function uses \code{Rcpp} for computational efficiency.
}
\seealso{
\code{\link{CGS_TwitterLDA}}, \code{\link{pred_TwitterLDA}}.
}
