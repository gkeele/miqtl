% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.DOQTL.to.HAPPY.R
\name{convert.full.DOQTL.array.to.HAPPY}
\alias{convert.full.DOQTL.array.to.HAPPY}
\title{Takes founder haplotype reconstructions as a 3D array (with the 36 probability states as one of the dimensions) 
and re-formats into a HAPPY-style genome cache}
\usage{
convert.full.DOQTL.array.to.HAPPY(
  full.array,
  map,
  map.locus_name.colname = "SNP_ID",
  map.chr.colname = "Chr",
  map.physical_dist.colname = "Mb_NCBI38",
  map.genetic_dist.colname = "cM",
  HAPPY.output.path,
  remove.chr.from.chr = FALSE,
  physical_dist.is.Mb = TRUE,
  allele.labels = LETTERS[1:8],
  chr = c(1:19, "X"),
  diplotype.order = c("DOQTL", "CC", "qtl2")
)
}
\arguments{
\item{full.array}{3D array that contains the founder probabilities. Should be dimension n x 36 x p, where n is the number of individuals
and p is the number of loci.}

\item{map}{The map file (which contains important information on the loci) loaded into R.}

\item{map.locus_name.colname}{DEFAULT: "SNP_ID". The column name in the map data that corresponds to locus/marker names.}

\item{map.chr.colname}{DEFAULT: "Chr". The column name in the map data that corresponds to chromosome.}

\item{map.physical_dist.colname}{DEFAULT: "Mb_NCBI38". The column name in the map data that corresponds to physical position.}

\item{map.genetic_dist.colname}{DEFAULT: "cM". The column name in the map data that corresponds to genetic position.}

\item{HAPPY.output.path}{The path to a directory that will be created as the HAPPY-format genome cache.}

\item{remove.chr.from.chr}{DEFAULT: FALSE. Option to remove "chr" from chromosome information. As in, change "chr1" to "1". The function
expects chromosome information to not include "chr" as a prefix.}

\item{physical_dist.is.Mb}{DEFAULT: TRUE. IF true, Mb will be converted to bp within the function.}

\item{allele.labels}{DEFAULT: NULL. Allows for specification of founder labels different from what is in the DO-QTL
output. The DEFAULT of NULL leads to using the labels from the DO-QTL output.}

\item{chr}{DEFAULT: c(1:19, "X"). Allows for specification of the chromosomes. DEFAULT is all the chromosomes from the mouse.}

\item{diplotype.order}{DEFAULT: "DOQTL". Specify the order of the diplotypes in the 3D array, so that the rotation to 
dosages is done correctly.}
}
\description{
This function produces a HAPPY-format genome cache from a founder haplotype reconstruction 3D array.
DO-QTL does not normally output this 3D array, but rather the dosages. However, this function requires
the full probabilities, not the dosages.
}
\examples{
convert.full.DOQTL.array.to.HAPPY()
}
