#' @export extract.mediation.qr
extract.mediation.qr <- function(genomecache, 
                                 id="SUBJECT.NAME", 
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
  pos <- as.numeric(gsub(pos, pattern="\\.[0-9]+\\.", replacement="", perl=TRUE))/1e6
  
  
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

#' @export mediation.scan.qr
mediation.scan.qr <- function(mediation.qr.object, 
                              data, 
                              phenotype,
                              return.allele.effects=FALSE,
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
  
  if(model == "full" & return.allele.effects){
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available with additive model\n",
        "Setting return.allele.effects to FALSE\n")
  }
  
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
    if(return.allele.effects){
      allele.effects[,i] <- get.allele.effects.from.fixef.eQTL(qr.alt=mediation.qr.object$qr.list[[mediator.index[i]]], 
                                                               y=y, 
                                                               founders=founders,
                                                               intercept.allele=mediation.qr.object$intercept.allele[i])
    }
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
                 allele.effects=allele.effects,
                 y=y,
                 formula=formula.string,
                 model.type=model,
                 p.value.method="ANOVA",
                 locus.effect.type="fixed",
                 n=length(y))
  return(output)
}

#' @export run.qr.permutation.threshold.mediation.scans
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
                                chr=chromatin.chr[chromatin.chr %in% chr.levels], 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=list(min=min.p,
                                               max=max.p))))
}

#' @export extract.chr.max.statistics.from.genomewide
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

#' @export extract.mediation.expression.qr
extract.mediation.expression.qr <- function(genomecache, 
                                            id="SUBJECT.NAME", 
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
  
  ## Extracting chromatin position
  pos <- gene.data[,2][gene.index]/1e6

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

