% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixef.eQTL.scan.R
\name{generate.qr.permutation.index.matrix}
\alias{generate.qr.permutation.index.matrix}
\title{Outputs a matrix of permuted indeces for permutation scans, allowing for the same permutations of individuals
across different outcomes.}
\usage{
generate.qr.permutation.index.matrix(
  qr.scan.object = NULL,
  n = NULL,
  num.samples,
  seed = 1
)
}
\arguments{
\item{qr.scan.object}{DEFAULT: NULL. Output object from scan.qr(). If NULL, function expects the number of individuals being
permutated.}

\item{n}{DEFAULT: NULL. Alternative to qr.scan.object. If NULL, function expects a qr.scan object.}

\item{num.samples}{The number of permutations of the index to create - ultimately the number of columns in the
output matrix.}

\item{seed}{DEFAULT: 1. Samplings of the index is a random process, thus a seed is necessary
to produce the same results over multiple runs and different machines.}
}
\description{
This function produces a matrix with columns that represent the permutations of the row index of the phenotype data.
This allows the same permutations of individuals to be used across different phenotypes. This approach is only statistically
valid when individuals are exchangeable.
}
\examples{
generate.qr.permutation.index.matrix()
}
