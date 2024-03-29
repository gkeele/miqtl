% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputed.snp.scan.h2lmm.R
\name{imputed.snp.scan.h2lmm}
\alias{imputed.snp.scan.h2lmm}
\title{Run a SNP-based genome scan from probabilities stored in a genome cache directory and founder strain
alleles stored .allele files}
\usage{
imputed.snp.scan.h2lmm(
  data,
  formula,
  K,
  allele.dir = NULL,
  genomecache,
  model = c("additive", "full"),
  use.par = "h2",
  chr = "all",
  brute = TRUE,
  use.fix.par = FALSE,
  just.these.loci = NULL,
  print.locus.fit = FALSE,
  use.progress.bar = TRUE,
  condition.loci = NULL,
  exclusion.freq = .Machine$double.eps,
  X.list = NULL,
  return.X.list = FALSE,
  ...
)
}
\arguments{
\item{data}{A data frame with outcome and potential covariates. Should also have individual IDs
that link to IDs in the genome cache  with a column named "SUBJECT.NAME".}

\item{formula}{An lm style formula with functions of outcome and covariates contained in data frame.}

\item{K}{DEFAULT: NULL. A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
based on SNP genotypes or the founder haplotype probabilities. Colnames and rownames should match
the SUBJECT.NAME column in the data frame. If no K matrix is specified, either lmer is used (if sparse random effects
are included in the formula) or a fixed effect model (equivalent to lm).}

\item{allele.dir}{DEFAULT: NULL. The path to the directory of .allele files that specify which SNP alleles correspond
to which founder haplotype. If NULL, X.list must be provided.}

\item{genomecache}{The path to the genome cache directory. The genome cache is a particularly structured
directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
store .RData files for the probabilities/dosages of each locus.}

\item{model}{DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
use of SNP dosages, and is most commonly used. The full option regresses the phenotype on the actual
genotype probabilities.}

\item{chr}{DEFAULT: "all". Specifies which chromosomes to scan.}

\item{brute}{DEFAULT: TRUE. During the optimization to find maximum likelihood parameters, this specifies checking the
boundaries of h2=0 and h2=1. Slightly less efficient, but otherwise the optimization procedure will not directly check
these values.}

\item{use.fix.par}{DEFAULT: TRUE. This specifies an approximate fitting of mixed effect model (Kang et al. 2009). Much
more efficient, as the optimization of h2 only needs to be performed once for the null model rather than every locus. 
Technically less powerful, though in practice it has proven to be almost equal to the exact procedure.}

\item{just.these.loci}{DEFAULT: NULL. Specifies a reduced set of loci to fit.}

\item{print.locus.fit}{DEFAULT: FALSE. If TRUE, prints out how many loci have been fit currently.}

\item{use.progress.bar}{DEFAULT: TRUE. If TRUE, a progress bar is used.}

\item{exclusion.freq}{DEFAULT: .Machine$double.eps. Loci with observed minor allele frequencies beneath
the specified value are removed from the scan.}

\item{X.list}{DEFAULT: NULL. This specifies the SNP-based design matrices for all the loci. If a scan
of the same population with the same markers has been performed, this option can save a lot of time.}

\item{return.X.list}{DEFAULT: FALSE. The scan procedure can return the list of design matrices for all
loci. It does increase the size of the output, though can then be used as an input for subsequent scans
of the same data.}
}
\description{
This function primarily takes a formula, data frame, genome cache, and directory of founder strain
.alleles files to run a genome scan.
}
\examples{
imputed.snp.scan.h2lmm()
}
