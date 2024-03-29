% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/condition_thru_residuals.R
\name{condition.out.locus.for.scan}
\alias{condition.out.locus.for.scan}
\title{Add the residuals to a data.frame after regressing out a specified QTL effect}
\usage{
condition.out.locus.for.scan(
  locus,
  new.outcome.name = "new.y",
  genomecache,
  data,
  formula,
  K = NULL,
  model = c("additive", "full"),
  locus.as.fixed = TRUE,
  use.par = "h2",
  use.multi.impute = TRUE,
  num.imp = 11,
  brute = TRUE,
  use.fix.par = TRUE,
  seed = 1,
  pheno.id = "SUBJECT.NAME",
  geno.id = "SUBJECT.NAME",
  weights = NULL,
  do.augment = FALSE,
  use.full.null = FALSE,
  added.data.points = 1,
  ...
)
}
\arguments{
\item{locus}{The locus for which the QTL effect will be regressed out. The locus must be contained in the
genomecache.}

\item{new.outcome.name}{DEFAULT: "new.y". The name of the new residual outcome variable to be included in the output data set.}

\item{genomecache}{The path to the genome cache directory. The genome cache is a particularly structured
directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
store .RData files for the probabilities/dosages of each locus.}

\item{data}{A data frame with outcome and potential covariates. Should also have IDs
that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
can be specified with pheno.id.}

\item{formula}{An lm style formula with functions of outcome and covariates contained in data frame.}

\item{K}{DEFAULT: NULL. A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
based on SNP genotypes or the founder haplotype probabilities. Colnames and rownames should match
the SUBJECT.NAME column in the data frame. If no K matrix is specified, either lmer is used (if sparse random effects
are included in the formula) or a fixed effect model (equivalent to lm).}

\item{model}{DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
diplotype probabilities.}

\item{locus.as.fixed}{DEFAULT: TRUE. If TRUE, the locus effect is fit as fixed effect. If FALSE, it is
fit as a random effect.}

\item{use.par}{DEFAULT: "h2". The parameterization of the likelihood to be used.}

\item{use.multi.impute}{DEFAULT: TRUE. This option specifies whether to use ROP or multiple imputations.}

\item{num.imp}{DEFAULT: 11. IF multiple imputations are used, this specifies the number of imputations to perform.}

\item{brute}{DEFAULT: TRUE. During the optimization to find maximum likelihood parameters, this specifies checking the
boundary of h2=0. Slightly less efficient, but otherwise the optimization procedure will not directly check
the boundary.}

\item{use.fix.par}{DEFAULT: TRUE. This specifies an approximate fitting of mixed effect model (Kang et al. 2009). Much
more efficient, as the optimization of h2 only needs to be performed once for the null model rather than every locus. 
Technically less powerful, though in practice it has proven to be almost equal to the exact procedure.}

\item{seed}{DEFAULT: 1. Multiple imputations involve a sampling process of the diplotypes, thus a seed is necessary
to produce the same results over multiple runs and different machines.}

\item{pheno.id}{DEFAULT: "SUBJECT.NAME". The is the individual-level ID that is associated with data points in the phenotype
data. Generally this should be unique for each data point.}

\item{geno.id}{DEFAULT: "SUBJECT.NAME". The default represents the situation that each genome is unique. Specifying some other
column allows for replicate genomes, such as in the CC or CC-RIX.}

\item{weights}{DEFAULT: NULL. If unspecified, individuals are equally weighted. This option allows for a weighted analysis 
when using the mean of multiple individuals with the same genome.}

\item{do.augment}{DEFAULT: FALSE. Augments the data with null observations for genotype groups. This is an approximately Bayesian 
approach to applying a prior to the data, and can help control highly influential data points.}

\item{use.full.null}{DEFAULT: FALSE. Draws augmented data points from the null model. This allows for the inclusion of null data points
that do not influence the estimation of other model parameters as much.}

\item{added.data.points}{DEFAULT: 1. If augment weights are being used, this specifies how many data points should be added in total.}
}
\description{
This function sets up conditional scans via regressing out QTL effects. The output is a new data.frame
with a residual outcome to be used as phenotype in subsequent genome scans.
}
\examples{
condition.out.locus.for.scan()
}
