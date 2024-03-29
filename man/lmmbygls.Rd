% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmmbygls.R
\name{lmmbygls}
\alias{lmmbygls}
\title{Run a haplotype-based genome scan from probabilities stored in a genome cache directory}
\usage{
lmmbygls(
  formula,
  data = NULL,
  y = NULL,
  X = NULL,
  K = NULL,
  eigen.K = NULL,
  fix.par = NULL,
  M = NULL,
  logDetV = NULL,
  weights = NULL,
  pheno.id = "SUBJECT.NAME",
  use.par = c("h2", "h2.REML"),
  brute = TRUE,
  subset,
  na.action,
  method = "qr",
  model = TRUE,
  contrasts = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{formula}{The lm/lmer-style formula for the model to fit. Variables must correspond to columns in the data.frame specified in the data
argument. If no data.frame is provided, and rather a y vector and X matrix, the formula is still passed to the resulting output fit object.}

\item{data}{DEFAULT: NULL. The data.frame that contains the variables included in the model. If no data.frame is specified, the expectation
is that a further process y vector and X matrix is provided. This allows for greater computational efficiency, side-stepping the need to pull
the quantities from the data.frame based on the formula.}

\item{y}{DEFAULT: NULL. y is an outcome vector that can be specified instead of a data.frame.}

\item{X}{DEFAULT: NULL. X is a design matrix that can be specified instead of a data.frame.}

\item{K}{DEFAULT: NULL. K is the covariance matrix, commonly a realized genetic relationship matrix or kinship matrix. NULL is interpreted as
no variance component should be fit. K should have row and column names that match the pheno.id column in the data.frame or the order of y and X.}

\item{eigen.K}{DEFAULT: NULL. Is the eigendecomposition of K. If NULL and K is non-null, it will be computed. Pre-computing in something like a
genome scan (as in scan.h2lmm()) can save computation time.}

\item{fix.par}{DEFAULT: NULL. This can be the ML or REML estimate of h2 from a null model. Using this within a genome scan results in an EMMAX-like
scan, which is the standard. It saves time because the parameter estimates will be analytical now rather than requiring an optimization step.}

\item{M}{DEFAULT: NULL. The GLS multiplier matrix. Useful in the context of genome scans with EMMAX-like specifications, as it saves computation time.}

\item{logDetV}{DEFAULT: NULL. Also useful within the context of a genome scan with EMMAX.}

\item{weights}{DEFAULT: NULL. Allows for the unstructured error to have diagonal but non-identity covariance structure. More intuitively, it allows 
observations to be differentially weighted. The default of NULL is equivalent to a vector of ones. Vector should have names that match pheno.id column
in data.frame or match the order of y and X.}

\item{pheno.id}{DEFAULT: "SUBJECT.NAME". The column in data that corresponds to the individual ID.}

\item{use.par}{DEFAULT: "h2". Specifies that the h2 should be optimized with respect to the likelihood (ML). Alternatively, "h2.REML" can be specified,
which optimizes h2 in terms of the residual likelihood (REML).}

\item{brute}{DEFAULT: TRUE. Internally optim() is used to find ML or REML estimates. It does not explicitly check the boundaries of h2. This option forces
it to check h2=0, which can occur.}

\item{subset}{Option to subset the data.}

\item{na.action}{Determines how NAs are handled.}

\item{method}{DEFAULT: "qr". Option used by lm.fit after GLS process is completed.}

\item{model}{DEFAULT: TRUE. Return model output. Can be turned off for efficiency within a scan.}

\item{contrasts}{DEFAULT: NULL.}

\item{verbose}{DEFAULT: FALSE.}
}
\description{
This function primarily takes a formula, data frame, and genome cache to run a genome scan.
}
\examples{
lmmbygls()
}
