#' Pre-compute QR decompositions from genome cache for efficient fixed effect model genome scan.
#'
#' This function calculates the QR decomposition for the null model, and all alternative models (per locus)
#' for a genome scan. This is only possible with a model that does not include any random effects.
#'
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param id DEFAULT: "SUBJECT.NAME". The is the individual-level ID that is associated with data points 
#' in the phenotype data. This should be unique for each data point.
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with pheno.id.
#' @param formula The right side of an lm-style formula with functions of covariates contained in data frame. First
#' symbol should be "~". Functions of the outcome will be specified in scan.qr().
#' @param model DEFAULT: "additive". Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities.
#' @param condition.loci DEFAULT: NULL. If loci are specified, they will be included in the null and alternative models.
#' The loci must be present in the genome cache. Alternatively, conditional scans can be accomplished by regressing out 
#' the loci effects with condition.out.locus.for.scan(), which does not require a new qr.object.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param just.these.loci DEFAULT: NULL. Specifies a reduced set of loci to fit. If loci is just one locus, the alternative model fit
#' will also be output as fit1.
#' @param use.progress.bar DEFAULT: TRUE. Results in a progress bar
#' @export
#' @examples extract.qr()
extract.qr <- function(genomecache, 
                       id="SUBJECT.NAME",
                       data, 
                       formula, 
                       model=c("additive", "full"), 
                       condition.loci=NULL,
                       chr="all", 
                       just.these.loci=NULL, 
                       use.progress.bar=TRUE){
  K <- NULL
  model <- model[1]
  
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  
  loci.chr <- h$getChromOfLocus(loci)
  if(chr[1] != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  subjects <- as.character(data[,id])
  X.0 <- model.matrix(formula, data=data)
  if(!is.null(condition.loci)){
    for(i in 1:length(condition.loci)){
      X.condition <- h$getLocusMatrix(condition.loci[i], model=model, subjects=subjects)
      keep.col <- 1:ncol(X.condition)
      max.column <- which.max(colSums(X.condition, na.rm=TRUE))[1]
      keep.col <- keep.col[keep.col != max.column]
      X.0 <- cbind(X.0, X.condition[,keep.col])
    }
  }
  qr.0 <- qr(X.0)
  
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(loci), style=3) }
  qr.list <- list()
  intercept.allele <- rep(NA, length(loci)) # For allele effects
  for(i in 1:length(loci)){
    X <- h$getLocusMatrix(loci[i], model=model, subjects=subjects)
    keep.col <- 1:ncol(X)
    max.column <- which.max(colSums(X, na.rm=TRUE))[1]
    intercept.allele[i] <- founders[max.column]
    keep.col <- keep.col[keep.col != max.column]
    X <- cbind(X.0, X[,keep.col])
    qr.list[[i]] <- qr(X)
    if(use.progress.bar){ setTxtProgressBar(pb, i) }
  }
  names(qr.list) <- loci
  
  qr.object <- list(qr.list=qr.list,
                    intercept.allele=intercept.allele,
                    condition.loci=condition.loci,
                    qr.0=qr.0,
                    chr=loci.chr,
                    pos=list(cM=h$getLocusStart(loci, scale="cM"),
                             Mb=h$getLocusStart(loci, scale="Mb")),
                    model=model,
                    founders=h$getFounders(),
                    subjects=subjects,
                    formula=Reduce(paste, deparse(formula)))
}

get.allele.effects.from.fixef.eQTL <- function(qr.alt, y, founders, intercept.allele, 
                                               center=TRUE, scale=FALSE){
  regression.effects <- qr.coef(qr=qr.alt, y=y)
  effects <- regression.effects[founders]
  intercept <- regression.effects["(Intercept)"]
  names(effects) <- founders
  
  effects <- effects + intercept
  effects[which(names(effects) == intercept.allele)] <- intercept
  return(as.vector(scale(effects, center=center, scale=scale)))
}

#' Runs genome scan from a QR decomposition object and phenotype data file.
#'
#' This function runs the genome scan based on a QR decomposition object and phenotype data file. This functionality only
#' works for fixed effect models.
#'
#' @param qr.object QR decomposition output from extract.qr().
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with id.
#' @param phenotype The column name (or function of column name) of a variable in data. This will become the outcome
#' of the genome scan.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param id DEFAULT: "SUBJECT.NAME". This is the individual-level ID that is associated with data points 
#' in the phenotype data. This should be unique for each data point.
#' @param just.these.loci DEFAULT: NULL. Specifies a reduced set of loci to fit. If loci is just one locus, the alternative model fit
#' will also be output as fit1.
#' @param debug.single.fit DEFAULT: FALSE. If TRUE, a browser() call is activated after the first locus is fit. This option
#' allows developers to more easily debug while still using the actual R package.
#' @param use.progress.bar DEFAULT: TRUE. Results in a progress bar while code runs.
#' @export
#' @examples scan.qr()
scan.qr <- function(qr.object, 
                    data, 
                    phenotype,
                    return.allele.effects=FALSE,
                    chr="all", 
                    id="SUBJECT.NAME",
                    just.these.loci=NULL,
                    debug.single.fit=FALSE, 
                    use.progress.bar=TRUE,
                    ...){
  model <- qr.object$model
  subjects <- qr.object$subjects
  founders <- qr.object$founders
  num.founders <- length(founders)
  loci <- names(qr.object$qr.list)
  loci.chr <- qr.object$chr
  rh.formula <- qr.object$formula

  if(model == "full" & return.allele.effects){
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available with additive model\n",
        "Setting return.allele.effects to FALSE\n")
  }
  
  ## Matching the subject order in the data with the qr object
  reorder <- match(subjects, data[,id])
  data <- data[reorder,]
  n <- nrow(data)

  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  formula.string <- paste(phenotype, rh.formula)
  formula <- formula(formula.string)

  allele.effects <- NULL
  p.vec <- rep(NA, length(loci))

  if(return.allele.effects){ 
    allele.effects <- matrix(NA, nrow=length(founders), ncol=length(loci),
                             dimnames=list(founders, loci))
  }
  y <- as.vector(model.frame(formula, data=data)[,1])
  names(y) <- subjects
  # Progress bar
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(loci), style=3) }
  for(i in 1:length(loci)){
    p.vec[i] <- get.f.stat.p.val(qr.alt=qr.object$qr.list[[i]], 
                                 qr.null=qr.object$qr.0, 
                                 y=y)
    if(is.nan(p.vec[i])){ p.vec[i] <- 1 }
    if(return.allele.effects){
      allele.effects[,i] <- get.allele.effects.from.fixef.eQTL(qr.alt=qr.object$qr.list[[i]], 
                                                               y=y, 
                                                               founders=founders,
                                                               intercept.allele=qr.object$intercept.allele[i])
    }
    if(debug.single.fit){ browser() }
    # Update progress bar
    if(use.progress.bar){ setTxtProgressBar(pb, i) }
  }
  names(p.vec) <- loci
  if(!is.null(qr.object$condition.loci)){
    formula.string <- paste(Reduce(paste, deparse(formula)), paste(qr.object$condition.loci, collapse="+"), sep="+")
  }
  output <- list(LOD=NULL,
                 p.value=p.vec,
                 df=NULL,
                 pos=list(Mb=qr.object$pos$Mb, 
                          cM=qr.object$pos$cM),
                 loci=loci, 
                 chr=loci.chr,
                 allele.effects=allele.effects,
                 y=y,
                 formula=formula.string,
                 model.type=model,
                 p.value.method="ANOVA",
                 locus.effect.type="fixed",
                 n=length(y))
  return(output)
}

#' Outputs a matrix of permuted indeces for permutation scans, allowing for the same permutations of individuals
#' across different outcomes.
#'
#' This function produces a matrix with columns that represent the permutations of the row index of the phenotype data.
#' This allows the same permutations of individuals to be used across different phenotypes. This approach is only statistically
#' valid when individuals are exchangeable.
#' 
#' @param qr.scan.object DEFAULT: NULL. Output object from scan.qr(). If NULL, function expects the number of individuals being
#' permutated.
#' @param n DEFAULT: NULL. Alternative to qr.scan.object. If NULL, function expects a qr.scan object.
#' @param num.samples The number of permutations of the index to create - ultimately the number of columns in the
#' output matrix.
#' @param seed DEFAULT: 1. Samplings of the index is a random process, thus a seed is necessary
#' to produce the same results over multiple runs and different machines.
#' @export
#' @examples generate.qr.permutation.index.matrix()
generate.qr.permutation.index.matrix <- function(qr.scan.object=NULL, 
                                                 n=NULL, 
                                                 num.samples, 
                                                 seed=1){
  if(!is.null(qr.scan.object)){
    n <- length(qr.scan.object$y)
  }
  
  set.seed(seed)
  perm.ind.matrix <- replicate(n=num.samples, sample(1:n, replace=FALSE))
  colnames(perm.ind.matrix) <- paste0("perm.", 1:num.samples)
  return(perm.ind.matrix)
}

#' Runs permutation scans for significance thresholds based on an index permutation matrix object, the phenotype data, and
#' pre-computed QR decompositions for all the models.
#'
#' This function runs permutation scans in order to calculate significance thresholds. Its results are most valid when
#' done in a fixed effect model setting and observations are exchangeable.
#' 
#' @param perm.ind.matrix Output object from generate.qr.permutation.index.matrix().
#' @param qr.object Output object from extract.qr().
#' @param phenotype The column name (or function of column name) of a variable in data. This will become the outcome
#' of the genome scan.
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with id.
#' @param keep.full.scans DEFAULT: FALSE. If TRUE, all p-values are kept from all loci. If FALSE, only the minimum p-value
#' is kept from each scan, greatly reducing size of output.
#' @param scan.index DEFAULT: NULL. If NULL, all permutations are run. Integer vector can be specified to run just a 
#' subset of the permutations.
#' @param id DEFAULT: "SUBJECT.NAME". This is the individual-level ID that is associated with data points 
#' in the phenotype data. This should be unique for each data point.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param just.these.loci DEFAULT: NULL. DEFAULT: NULL. Specifies a reduced set of loci to fit. If loci is just one locus, the alternative model fit
#' will also be output as fit1.
#' @param use.progress.bar DEFAULT: FALSE. Results in a progress bar while code runs.
#' @param report.perm.scans DEFAULT: FALSE. Outputs the completion of permutation scans.
#' @export
#' @examples run.qr.permutation.threshold.scans()
run.qr.permutation.threshold.scans <- function(perm.ind.matrix, 
                                               qr.object,
                                               phenotype, 
                                               data,
                                               keep.full.scans=FALSE, 
                                               scan.index=NULL, 
                                               id="SUBJECT.NAME",
                                               chr="all", 
                                               just.these.loci=NULL, 
                                               use.progress.bar=FALSE,
                                               report.perm.scans=FALSE,
                                               ...){
  if(is.null(scan.index)){ scan.index <- 1:ncol(perm.ind.matrix) }
  
  loci <- names(qr.object$qr.list)
  rh.formula <- qr.object$formula
  loci.chr <- qr.object$chr
  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  formula.string <- paste(phenotype, rh.formula)
  formula <- formula(formula.string)
  
  full.p <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- matrix(NA, nrow=length(scan.index), ncol=length(loci))
    colnames(full.p) <- loci
    these.pos <- list(Mb=qr.object$pos$Mb[loci],
                      cM=qr.object$pos$cM[loci])
  }
  min.p <- rep(NA, length(scan.index))
  for(i in 1:length(scan.index)){
    data$new_y <- model.frame(formula, data=data)[perm.ind.matrix[,scan.index[i]], 1]

    this.scan <- scan.qr(qr.object=qr.object, data=data, 
                         phenotype="new_y", id=id, chr=chr, 
                         return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                         ...)
    if(keep.full.scans){
      full.p[i,] <- this.scan$p.value
    }
    min.p[i] <-  min(this.scan$p.value)
    if (report.perm.scans) {
      cat("\n", "Threshold scan: index", scan.index[i], "complete ---------- final index of this run:", scan.index[length(scan.index)], "\n")
    }
  }
  return(list(full.results=list(LOD=NULL,
                                p.value=full.p,
                                chr=loci.chr, 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=min.p)))
}

