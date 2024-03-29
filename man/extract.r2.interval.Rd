% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp_r2.R
\name{extract.r2.interval}
\alias{extract.r2.interval}
\title{Calculate the r^2 (squared correlation coefficient) between the genotype at all loci on a chromosome 
and a specified SNP marker.}
\usage{
\method{extract}{r2.interval}(scan.object, r2.scan.object, r2.level = 0.6)
}
\arguments{
\item{scan.object}{An SNP scan object produced by imputed.snp.scan.h2lmm().}

\item{r2.scan.object}{An r^2 object produced by pairwise.cor.snp.scan().}

\item{r2.level}{DEFAULT: 0.6. The r^2 cutpoint. Returns the position interval includes loci with
r^2 \eqn{\ge} r2.level.}
}
\description{
This function primarily takes a formula, data frame, genome cache, and directory of founder strain
.alleles files to calculate the pairwise r^2 between all individual markers and a specified
locus, likely the peak SNP.
}
\examples{
extract.r2.interval()
}
