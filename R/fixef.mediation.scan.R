#' Pre-compute QR decompositions from genome cache for efficient fixed effect model mediation scans.
#'
#' This function calculates the QR decomposition for all null model and alternative models (per chromatin mediators)
#' for a genome chromatin mediation scan. This is only possible with a model that does not include any random effects.
#' Currently designed for chromatin mediation (measured by ATAC-seq), though could be generalized.
#'
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param id DEFAULT: "SUBJECT.NAME". The is the individual-level ID that is associated with data points 
#' in the phenotype data. This should be unique for each data point.
#' @param pos.is.bp DEFAULT: TRUE. If the pos variable in data is bp, use the default. If it is Mb, then
#' set to FALSE.
#' @param locus The marker or interval at which the QTL is detected. For the mediation scan, the detected
#' QTL is fit in the alternative model for each mediator.
#' @param data A data frame with mediators and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with pheno.id. Currently designed to filter data to the specified covariates and the 
#' chromatin variables.
#' @param formula The right side of an lm-style formula with functions of covariates contained in data frame. First
#' symbol should be "~". Functions of the outcome will be specified in mediation.scan.qr().
#' @param model DEFAULT: "additive". Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities.
#' @param condition.loci DEFAULT: NULL. If loci are specified, they will be included in the null and alternative models.
#' The loci must be present in the genome cache. Alternatively, conditional scans can be accomplished by regressing out 
#' the loci effects with condition.out.locus.for.scan(), which does not require a new qr.object.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param use.progress.bar DEFAULT: TRUE. Results in a progress bar being displayed.
#' @export extract.mediation.qr
#' @examples extract.mediation.qr()
extract.mediation.qr <- function(genomecache, 
                                 id="SUBJECT.NAME", 
                                 pos.is.bp = TRUE,
                                 locus,
                                 data, 
                                 formula, 
                                 model=c("additive", "full"), 
                                 condition.loci=NULL,
                                 chr="all", 
                                 use.progress.bar=TRUE){
  K <- NULL
  model <- model[1]
  
  ## Extracting chromatin phenotypes
  chromatin <- grep(names(data), pattern="^chr", perl=TRUE, value=TRUE)
  ## Extracting chr for chromatin phenotypes
  chromatin.chr <- gsub(chromatin, pattern="\\.[0-9]+\\.[0-9]+$", replacement="", perl=TRUE)
  chromatin.chr <- gsub(chromatin.chr, pattern="^chr", replacement="", perl=TRUE)
  ## Extracting chromatin position
  pos <- gsub(chromatin, pattern="^chr[0-9X]+\\.", replacement="", perl=TRUE)
  if (pos.is.bp) {
    pos <- as.numeric(gsub(pos, pattern="\\.[0-9]+\\.", replacement="", perl=TRUE))/1e6
  }
  
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  
  locus.chr <- h$getChromOfLocus(locus)
  if(chr[1] != "all"){
    chromatin <- chromatin[chromatin.chr %in% chr]
    chromatin.chr <- chromatin.chr[chromatin.chr %in% chr]
  }
  
  subjects <- as.character(data[,id])
  X.0 <- model.matrix(formula, data=data)
  X.locus <- h$getLocusMatrix(locus, model=model, subjects=subjects)
  keep.col <- 1:ncol(X.locus)
  max.column <- which.max(colSums(X.locus, na.rm=TRUE))[1]
  keep.col <- keep.col[keep.col != max.column]
  if(!is.null(condition.loci)){
    for(i in 1:length(condition.loci)){
      X.condition <- h$getLocusMatrix(condition.loci[i], model=model, subjects=subjects)
      keep.col <- 1:ncol(X.condition)
      max.column <- which.max(colSums(X.condition, na.rm=TRUE))[1]
      keep.col <- keep.col[keep.col != max.column]
      X.0 <- cbind(X.0, X.condition[,keep.col])
    }
  }
  
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(chromatin), style=3) }
  qr.0.list <- qr.alt.list <- list()
  intercept.allele <- rep(founders[max.column], length(chromatin)) # For allele effects
  for(i in 1:length(chromatin)){
    X <- rint(data[,chromatin[i]])
    
    this.X.0 <- cbind(X.0, X)
    this.X <- cbind(X.0, X.locus, X)
    qr.0.list[[i]] <- qr(this.X.0)
    qr.alt.list[[i]] <- qr(this.X)
    if(use.progress.bar){ setTxtProgressBar(pb, i) }
  }
  names(qr.alt.list) <- names(qr.0.list) <- chromatin
  
  qr.object <- list(qr.alt.list=qr.alt.list,
                    qr.0.list=qr.0.list,
                    X.0=X.0,
                    X.locus=X.locus,
                    intercept.allele=intercept.allele,
                    condition.loci=condition.loci,
                    chr=chromatin.chr,
                    locus=locus,
                    locus.chr=locus.chr,
                    pos=list(cM=NA,
                             Mb=pos),
                    model=model,
                    founders=founders,
                    subjects=subjects,
                    formula=Reduce(paste, deparse(formula)))
}

#' Runs genome-wide mediation scan from a QR decomposition object based on a specified QTL and phenotype data file.
#'
#' This function runs the genome scan based on a QR decomposition object and phenotype data file. This functionality only
#' works for fixed effect models. Currently written for chromatin mediation of local-eQTL, though could easily be 
#' generalized to any pairing of mediator and outcome.
#'
#' @param mediation.qr.object Output object from extract.mediation.qr()
#' @param data A data frame with outcome. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with pheno.id.
#' @param phenotype The column name (or function of column name) of a variable in data. This will become the outcome
#' of the mediation scan.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param id DEFAULT: "SUBJECT.NAME". The is the individual-level ID that is associated with data points 
#' in the phenotype data. This should be unique for each data point.
#' @param debug.single.fit DEFAULT: FALSE. If TRUE, a browser() call is activated after the first locus is fit. This option
#' allows developers to more easily debug while still using the actual R package.
#' @param use.progress.bar DEFAULT: TRUE. Results in a progress bar being displayed.
#' @export mediation.scan.qr
#' @examples mediation.scan.qr()
mediation.scan.qr <- function(mediation.qr.object, 
                              data, 
                              phenotype,
                              chr="all", 
                              id="SUBJECT.NAME",
                              debug.single.fit=FALSE, 
                              use.progress.bar=TRUE,
                              ...){
  model <- mediation.qr.object$model
  subjects <- mediation.qr.object$subjects
  founders <- mediation.qr.object$founders
  num.founders <- length(founders)
  mediator <- names(mediation.qr.object$qr.alt.list)
  mediator.chr <- mediation.qr.object$chr
  mediator.pos <- mediation.qr.object$pos$Mb
  rh.formula <- mediation.qr.object$formula
  
  ## Matching the subject order in the data with the qr object
  reorder <- match(subjects, data[,id])
  data <- data[reorder,]
  n <- nrow(data)
  
  if (all.vars(formula(paste0(phenotype, "~1")))[1] %in% mediator) {
   remove.outcome <- which(all.vars(formula(paste0(phenotype, "~1")))[1] == mediator)
   mediator <- mediator[-remove.outcome]
   mediator.chr <- mediator.chr[-remove.outcome]
   mediator.pos <- mediator.pos[-remove.outcome]
  }
  if (chr[1] != "all") {
    mediator <- mediator[mediator.chr %in% chr]
    mediator.chr <- mediator.chr[mediator.chr %in% chr]
    mediator.pos <- mediator.pos[mediator.chr %in% chr]
  }
  
  formula.string <- paste(phenotype, rh.formula)
  formula <- formula(formula.string)
  
  allele.effects <- NULL
  p.vec <- rep(NA, length(mediator))
  
  if(return.allele.effects){ 
    allele.effects <- matrix(NA, nrow=length(founders), ncol=length(mediator),
                             dimnames=list(founders, mediator))
  }
  y <- as.vector(model.frame(formula, data=data)[,1])
  names(y) <- subjects
  # Progress bar
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(mediator), style=3) }
  mediator.index <- unlist(sapply(1:length(mediator), function(x) which(mediator[x] == names(mediation.qr.object$qr.alt.list))))
    
  for(i in 1:length(mediator)){
    p.vec[i] <- get.f.stat.p.val(qr.alt=mediation.qr.object$qr.alt.list[[mediator.index[i]]], 
                                 qr.null=mediation.qr.object$qr.0.list[[mediator.index[i]]], 
                                 y=y)
    if(is.nan(p.vec[i])){ p.vec[i] <- 1 }
    if(debug.single.fit){ browser() }
    # Update progress bar
    if(use.progress.bar){ setTxtProgressBar(pb, i) }
  }
  names(p.vec) <- mediator
  if(!is.null(mediation.qr.object$condition.loci)){
    formula.string <- paste(Reduce(paste, deparse(formula)), paste(mediation.qr.object$condition.loci, collapse="+"), sep="+")
  }
  output <- list(LOD=NULL,
                 p.value=p.vec,
                 df=NULL,
                 pos=list(Mb=mediator.pos, 
                          cM=mediation.qr.object$pos$cM),
                 loci=mediator, 
                 chr=mediator.chr,
                 y=y,
                 formula=formula.string,
                 model.type=model,
                 p.value.method="ANOVA",
                 locus.effect.type="fixed",
                 locus=mediation.qr.object$locus,
                 n=length(y))
  return(output)
}

#' Runs mediation permutation scans for significance thresholds based on an index permutation matrix object, the phenotype data, and
#' pre-computed QR decompositions for all the models.
#'
#' This function runs permutation scans in order to calculate significance thresholds for mediation analysis. Its results are most valid when
#' done in a fixed effect model setting and observations are exchangeable.
#' 
#' @param perm.ind.matrix Output object from generate.qr.permutation.index.matrix().
#' @param mediation.qr.object Output object from extract.mediation.qr().
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
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
#' @param use.progress.bar DEFAULT: FALSE. Results in a progress bar while code runs.
#' @param pos.is.bp DEFAULT: TRUE. If the pos variable in data is bp, use the default. If it is Mb, then
#' set to FALSE.
#' @export run.qr.permutation.threshold.mediation.scans
#' @examples run.qr.permutation.threshold.mediation.scans()
run.qr.permutation.threshold.mediation.scans <- function(perm.ind.matrix, 
                                                         mediation.qr.object, 
                                                         genomecache,
                                                         phenotype, 
                                                         data,
                                                         keep.full.scans=FALSE, 
                                                         scan.index=NULL, 
                                                         id="SUBJECT.NAME",
                                                         chr="all", 
                                                         use.progress.bar=TRUE,
                                                         pos.is.bp = TRUE,
                                                         ...){
  
  if(is.null(scan.index)){ scan.index <- 1:ncol(perm.ind.matrix) }
  
  model <- mediation.qr.object$model
  condition.loci <- mediation.qr.object$condition.loci
  locus <- mediation.qr.object$locus
  
  chromatin <- names(mediation.qr.object$qr.0.list)
  rh.formula <- mediation.qr.object$formula
  chromatin.chr <- mediation.qr.object$chr
  if(chr[1] != "all"){
    chromatin <- chromatin[chromatin.chr %in% chr]
  }
  formula.string <- paste(phenotype, rh.formula, "+", id)
  formula <- formula(formula.string)
  
  if(chr == "all"){ chr.levels <- unique(chromatin.chr) }
  else{ chr.levels <- chr }
  
  full.p <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- matrix(NA, nrow=length(scan.index), ncol=length(chromatin))
    colnames(full.p) <- chromatin
    these.pos <- list(Mb=mediation.qr.object$pos$Mb[chromatin.chr %in% chr.levels],
                      cM=mediation.qr.object$pos$cM[chromatin.chr %in% chr.levels])
  }
  min.p <- max.p <- rep(NA, length(scan.index))
  y <- model.frame(formula, data=data)[,1]
  data <- data[,!grepl(pattern="^gene_", perl=TRUE, x=colnames(data))]
  permute.var <- !(colnames(data) %in% all.vars(formula)[-1])
  for(i in 1:length(scan.index)){
    ## Permuting all variables but covariates
    this.data <- data.frame(y=y, data[perm.ind.matrix[,scan.index[i]], permute.var], data[, !permute.var])
    this.mediation.qr.object <- extract.mediation.qr(genomecache=genomecache, id=id,
                                                     data=this.data, formula=as.formula(rh.formula), 
                                                     model=model, condition.loci=condition.loci,
                                                     chr=chr, locus=locus, use.progress.bar=FALSE,
                                                     pos.is.bp = pos.is.bp)
    this.scan <- mediation.scan.qr(mediation.qr.object=this.mediation.qr.object, data=this.data, 
                                   phenotype="y", id=id, chr=chr, 
                                   return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                                   ...)
    if(keep.full.scans){
      full.p[i,] <- this.scan$p.value
    }
    min.p[i] <-  min(this.scan$p.value)
    max.p[i] <- max(this.scan$p.value)
    cat("\n", "Threshold scan: index", scan.index[i], "complete ---------- final index of this run:", scan.index[length(scan.index)], "\n")
  }
  return(list(full.results=list(LOD=NULL,
                                p.value=full.p,
                                chr=chromatin.chr[chromatin.chr %in% chr.levels], 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=list(min=min.p,
                                               max=max.p))))
}

#' Extracts chromosome-level statistics from the full output of genome-wide permutation scans, allowing chromosome-wide and genome-wide significance
#' to be determined from the same output.
#'
#' This function extracts the associations for a specified chromosome from genome-wide permutation scans. This allows both genome-wide and chromosome-wide
#' thresholds to be determined from the same output.
#' 
#' @param full.perm.scans Output object from generate.qr.permutation.index.matrix().
#' @param chr The chromosome that will have its association scores extracted from the full.perm.scans, essentially
#' creating permutation scans specific to the chromosome.
#' @param use.lod DEFAULT: FALSE. If TRUE, LOD scores are recorded. If FALSE, p-values are recorded, either based on the likelihood ratio test or ANOVA.
#' @export extract.chr.max.statistics.from.genomewide
#' @examples extract.chr.max.statistics.from.genomewide()
extract.chr.max.statistics.from.genomewide <- function(full.perm.scans, 
                                                       chr, 
                                                       use.lod=FALSE){
  statistic.type <- ifelse(use.lod, "LOD", "p.value")
  full.results <- full.perm.scans$full.results[[statistic.type]][,full.perm.scans$full.results$chr == chr]
  
  return(list(full.results=NULL,
              max.statistics=list(LOD=NULL,
                                  p.value=list(min=apply(full.results, 1, function(x) min(x, na.rm=TRUE)),
                                               max=apply(full.results, 1, function(x) max(x, na.rm=TRUE))))))
}

#' Pre-compute QR decompositions from genome cache for efficient fixed effect model mediation scans.
#'
#' This function calculates the QR decomposition for all null model and alternative models (per gene expression mediators)
#' for a genome expression mediation scan. This is only possible with a model that does not include any random effects.
#' Currently designed for gene expression mediation (measured by RNA-seq), though could be generalized, essentially
#' consolidating extract.mediation.expression.qr() and extract.mediation.qr().
#'
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param id DEFAULT: "SUBJECT.NAME". The is the individual-level ID that is associated with data points 
#' in the phenotype data. This should be unique for each data point.
#' @param pos.is.bp DEFAULT: TRUE. If the pos variable in gene.data is bp, use the default. If it is Mb, then
#' set to FALSE.
#' @param locus The marker or interval at which the QTL is detected. For the mediation scan, the detected
#' QTL is fit in the alternative model for each mediator.
#' @param data A data frame with mediators and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with pheno.id. Currently designed to filter data to the specified covariates and the 
#' gene expression variables.
#' @param gene.data Gene annotation data for the gene expression variables in data. Specifically, the gene transcription
#' start site is needed, named gene.start in gene.data.
#' @param formula The right side of an lm-style formula with functions of covariates contained in data frame. First
#' symbol should be "~". Functions of the outcome will be specified in mediation.scan.qr().
#' @param model DEFAULT: "additive". Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities.
#' @param condition.loci DEFAULT: NULL. If loci are specified, they will be included in the null and alternative models.
#' The loci must be present in the genome cache. Alternatively, conditional scans can be accomplished by regressing out 
#' the loci effects with condition.out.locus.for.scan(), which does not require a new qr.object.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param use.progress.bar DEFAULT: TRUE. Results in a progress bar being displayed.
#' @export extract.mediation.expression.qr
#' @examples extract.mediation.expression.qr()
extract.mediation.expression.qr <- function(genomecache, 
                                            id="SUBJECT.NAME", 
                                            pos.is.bp = TRUE,
                                            locus,
                                            data, 
                                            gene.data,
                                            formula, 
                                            model=c("additive", "full"), 
                                            condition.loci=NULL,
                                            chr="all", 
                                            use.progress.bar=TRUE){
  K <- NULL
  model <- model[1]
  
  
  ## Extracting gene expressions
  genes <- grep(names(data), pattern="^gene", perl=TRUE, value=TRUE)
  ## Extracting chr for expression phenotypes
  ref.genes <- paste0("gene_", gene.data[,4])
  ref.genes <- chartr(old="()-", new="...", ref.genes)
  
  gene.index <- unlist(sapply(1:length(genes), function(x) which(genes[x] == ref.genes)))
  gene.chr <- gene.data[,1][gene.index]
  gene.chr <- gsub(gene.chr, pattern="^chr", replacement="", perl=TRUE)
  
  ## Extracting gene position
  pos <- gene.data[,2][gene.index]
  if (pos.is.bp) {
    pos <- pos/1e6
  }

  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  
  locus.chr <- h$getChromOfLocus(locus)
  if(chr[1] != "all"){
    genes <- genes[gene.chr %in% chr]
    gene.chr <- gene.chr[gene.chr %in% chr]
    pos <- pos[gene.chr %in% chr]
  }
  
  subjects <- as.character(data[,id])
  X.0 <- model.matrix(formula, data=data)
  X.locus <- h$getLocusMatrix(locus, model=model, subjects=subjects)
  keep.col <- 1:ncol(X.locus)
  max.column <- which.max(colSums(X.locus, na.rm=TRUE))[1]
  keep.col <- keep.col[keep.col != max.column]
  if(!is.null(condition.loci)){
    for(i in 1:length(condition.loci)){
      X.condition <- h$getLocusMatrix(condition.loci[i], model=model, subjects=subjects)
      keep.col <- 1:ncol(X.condition)
      max.column <- which.max(colSums(X.condition, na.rm=TRUE))[1]
      keep.col <- keep.col[keep.col != max.column]
      X.0 <- cbind(X.0, X.condition[,keep.col])
    }
  }
  
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(genes), style=3) }
  qr.0.list <- qr.alt.list <- list()
  intercept.allele <- rep(founders[max.column], length(genes)) # For allele effects
  for(i in 1:length(genes)){
    X <- rint(data[,genes[i]])
    
    this.X.0 <- cbind(X.0, X)
    this.X <- cbind(X.0, X.locus, X)
    qr.0.list[[i]] <- qr(this.X.0)
    qr.alt.list[[i]] <- qr(this.X)
    if(use.progress.bar){ setTxtProgressBar(pb, i) }
  }
  names(qr.alt.list) <- names(qr.0.list) <- genes
  
  qr.object <- list(qr.alt.list=qr.alt.list,
                    qr.0.list=qr.0.list,
                    X.0=X.0,
                    X.locus=X.locus,
                    intercept.allele=intercept.allele,
                    condition.loci=condition.loci,
                    chr=gene.chr,
                    locus=locus,
                    locus.chr=locus.chr,
                    pos=list(cM=NA,
                             Mb=pos),
                    model=model,
                    founders=founders,
                    subjects=subjects,
                    formula=Reduce(paste, deparse(formula)))
}

#' Runs mediation permutation scans for significance thresholds based on an index permutation matrix object, the phenotype data, and
#' pre-computed QR decompositions for all the models. Currently designed for mediation through gene expression.
#'
#' This function runs permutation scans in order to calculate significance thresholds for mediation analysis. Its results are most valid when
#' done in a fixed effect model setting and observations are exchangeable.
#' 
#' @param perm.ind.matrix Output object from generate.qr.permutation.index.matrix().
#' @param mediation.qr.object Output object from extract.mediation.expression.qr().
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
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
#' @param use.progress.bar DEFAULT: FALSE. Results in a progress bar while code runs.
#' @param pos.is.bp DEFAULT: TRUE. If the pos variable in data is bp, use the default. If it is Mb, then
#' set to FALSE.
#' @export run.qr.permutation.threshold.mediation.expression.scans
#' @examples run.qr.permutation.threshold.mediation.expression.scans()
run.qr.permutation.threshold.mediation.expression.scans <- function(perm.ind.matrix, 
                                                                    mediation.qr.object, 
                                                                    genomecache,
                                                                    phenotype, 
                                                                    data,
                                                                    gene.data,
                                                                    keep.full.scans=FALSE, 
                                                                    scan.index=NULL, 
                                                                    id="SUBJECT.NAME",
                                                                    chr="all", 
                                                                    pos.is.bp=TRUE,
                                                                    use.progress.bar=TRUE,
                                                                    ...){
  
  if(is.null(scan.index)){ scan.index <- 1:ncol(perm.ind.matrix) }
  
  model <- mediation.qr.object$model
  condition.loci <- mediation.qr.object$condition.loci
  locus <- mediation.qr.object$locus
  
  genes <- names(mediation.qr.object$qr.0.list)
  rh.formula <- mediation.qr.object$formula
  genes.chr <- mediation.qr.object$chr
  if(chr[1] != "all"){
    genes <- genes[genes.chr %in% chr]
  }
  formula.string <- paste(phenotype, rh.formula, "+", id)
  formula <- formula(formula.string)
  
  if(chr == "all"){ chr.levels <- unique(genes.chr) }
  else{ chr.levels <- chr }
  
  full.p <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- matrix(NA, nrow=length(scan.index), ncol=length(genes))
    colnames(full.p) <- genes
    these.pos <- list(Mb=mediation.qr.object$pos$Mb[genes.chr %in% chr.levels],
                      cM=mediation.qr.object$pos$cM[genes.chr %in% chr.levels])
  }
  min.p <- max.p <- rep(NA, length(scan.index))
  y <- model.frame(formula, data=data)[,1]
  permute.var <- !(colnames(data) %in% all.vars(formula)[-1])
  for(i in 1:length(scan.index)){
    ## Permuting all variables but covariates
    this.data <- data.frame(y=y, data[perm.ind.matrix[,scan.index[i]], permute.var], data[, !permute.var])
    this.mediation.qr.object <- extract.mediation.expression.qr(genomecache=genomecache, id=id,
                                                                data=this.data, 
                                                                gene.data=gene.data, 
                                                                formula=as.formula(rh.formula),
                                                                pos.is.bp=pos.is.bp,
                                                                model=model, condition.loci=condition.loci,
                                                                chr=chr, locus=locus, use.progress.bar=FALSE)
    this.scan <- mediation.scan.qr(mediation.qr.object=this.mediation.qr.object, data=this.data, 
                                   phenotype="y", id=id, chr=chr, 
                                   return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                                   ...)
    if(keep.full.scans){
      full.p[i,] <- this.scan$p.value
    }
    min.p[i] <-  min(this.scan$p.value)
    max.p[i] <- max(this.scan$p.value)
    cat("\n", "Threshold scan: index", scan.index[i], "complete ---------- final index of this run:", scan.index[length(scan.index)], "\n")
  }
  return(list(full.results=list(LOD=NULL,
                                p.value=full.p,
                                chr=genes.chr[genes.chr %in% chr.levels], 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=list(min=min.p,
                                               max=max.p))))
}


