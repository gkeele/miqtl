#' Calculate the r^2 (squared correlation coefficient) between the genotype at all loci on a chromosome 
#' and a specified SNP marker.
#'
#' This function primarily takes a formula, data frame, genome cache, and directory of founder strain
#' .alleles files to calculate the pairwise r^2 between all individual markers and a specified
#' locus, likely the peak SNP.
#'
#' @param data A data frame with outcome and potential covariates. Should also have individual IDs
#' that link to IDs in the genome cache  with a column named "SUBJECT.NAME".
#' @param formula An lm style formula with functions of outcome and covariates contained in data frame.
#' @param K DEFAULT: NULL. A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
#' based on SNP genotypes or the founder haplotype probabilities. Colnames and rownames should match
#' the SUBJECT.NAME column in the data frame. If no K matrix is specified, either lmer is used (if sparse random effects
#' are included in the formula) or a fixed effect model (equivalent to lm).
#' @param allele.dir The path to the directory of .allele files that specify which SNP alleles correspond
#' to which founder haplotype. .allele files for a format used by HAPPY.
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param model DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of SNP dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' genotype probabilities.
#' @param chr Specifies which individual chromosomes to scan.
#' @param point.locus The locus to calculate all pairwise r^2 between other loci on the chromosome. Often
#' the peak SNP from a scan.
#' @param just.these.loci DEFAULT: NULL. Specifies a reduced set of loci to fit.
#' @param print.progress DEFAULT: FALSE. If TRUE, prints out how many loci have been fit currently.
#' @param exclusion.freq DEFAULT: .Machine$double.eps. Loci with observed minor allele frequencies beneath
#' the specified value are removed from the scan.
#' @param X.list DEFAULT: NULL. This specifies the SNP-based design matrices for all the loci. If a scan
#' of the same population with the same markers has been performed, this option can save a lot of time.
#' @export
#' @examples pairwise.cor.snp.scan()
pairwise.cor.snp.scan <- function(data, formula, K,
                                  allele.dir, genomecache,
                                  model=c("additive", "full"), chr, point.locus,
                                  just.these.loci=NULL,
                                  print.progress=FALSE,
                                  exclusion.freq=.Machine$double.eps,
                                  X.list,
                                  ...){
  model <- model[1]
  
  do.augment <- FALSE # Necessary to work with old support_functions.R
  
  mapping.matrix <- straineff.mapping.matrix()
  
  these.chr <- chr
  # Genomecache for imputation
  h <- DiploprobReader$new(genomecache)
  
  loci <- names(X.list)
  
  founders <- h$getFounders()
  cache.subjects <- rownames(h$getLocusMatrix(locus=loci[1], model="additive"))
  loci.chr <- h$getChromOfLocus(loci)
  
  loci <- loci[loci.chr == these.chr]
  loci.chr <- loci.chr[loci.chr == these.chr]
  
  data.and.K <- make.processed.data(formula=formula, data=data, K=K, cache.subjects=cache.subjects, impute.on="SUBJECT.NAME")
  data <- data.and.K$data
  
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  null.data <- data
  null.K <- K
  
  # Grabbing the design matrix (vector) of the point SNP
  point.X <- X.list[[point.locus]][as.character(null.data$SUBJECT.NAME),, drop=FALSE]
  
  r2.vec <- rep(0, length(loci))
  for(i in 1:length(loci)){
    X <- X.list[[loci[i]]][as.character(null.data$SUBJECT.NAME),, drop=FALSE]
    
    r2.vec[i] <- cor(point.X, X)^2
    
    if(print.progress){ cat(paste("locus", i, "out of", length(loci)), "\n") }
  }
  output <- list(r2=r2.vec,
                 point.locus=point.locus,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=chr,
                 model.type=model)
  return(output)
}

#' Calculate the r^2 (squared correlation coefficient) between the genotype at all loci on a chromosome 
#' and a specified SNP marker.
#'
#' This function primarily takes a formula, data frame, genome cache, and directory of founder strain
#' .alleles files to calculate the pairwise r^2 between all individual markers and a specified
#' locus, likely the peak SNP.
#'
#' @param scan.object An SNP scan object produced by imputed.snp.scan.h2lmm().
#' @param r2.scan.object An r^2 object produced by pairwise.cor.snp.scan().
#' @param r2.level DEFAULT: 0.6. The r^2 cutpoint. Returns the position interval includes loci with
#'r^2 \eqn{\ge} r2.level.
#' @export
#' @examples extract.r2.interval()
extract.r2.interval <- function(scan.object, r2.scan.object, r2.level=0.6){
  high.r2 <- r2.scan.object$r2[which(r2.scan.object$r2 > r2.level)]
  high.r2.loci <- r2.scan.object$loci[which(r2.scan.object$r2 > r2.level)]
  high.r2.pvalue <- scan.object$p.value[which(scan.object$loci %in% high.r2.loci)]
  high.r2.pos.Mb <- scan.object$pos$Mb[which(scan.object$loci %in% high.r2.loci)]
  high.r2.pos.cM <- scan.object$pos$cM[which(scan.object$loci %in% high.r2.loci)]
  max.pos.index <- which.max(high.r2.pos.cM)
  min.pos.index <- which.min(high.r2.pos.cM)
  ## Interval
  interval.width.cM <- high.r2.pos.cM[max.pos.index] - high.r2.pos.cM[min.pos.index]
  interval.width.Mb <- high.r2.pos.Mb[max.pos.index] - high.r2.pos.Mb[min.pos.index]
  ## Boundary markers
  ub.marker <- high.r2.loci[max.pos.index]
  lb.marker <- high.r2.loci[min.pos.index]
  ub.cM <- high.r2.pos.cM[max.pos.index]
  lb.cM <- high.r2.pos.cM[min.pos.index]
  ub.Mb <- high.r2.pos.Mb[max.pos.index]
  lb.Mb <- high.r2.pos.Mb[min.pos.index]
  
  results <- data.frame(peak.marker=r2.scan.object$point.locus, r2.level=r2.level,
                        ub.marker=ub.marker, ub.cM=round(ub.cM, 2), ub.Mb=round(ub.Mb, 2),
                        lb.marker=lb.marker, lb.cM=round(lb.cM, 2), lb.Mb=round(lb.Mb, 2),
                        interval.width.cM=round(interval.width.cM, 2), interval.width.Mb=round(interval.width.Mb, 2))
  return(results)
}
