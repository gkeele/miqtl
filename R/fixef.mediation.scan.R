#' @export
extract.mediation.qr <- function(genomecache, id="SUBJECT.NAME", locus,
                                 data, formula, model=c("additive", "full"), condition.loci=NULL,
                                 chr="all", use.progress.bar=TRUE){
  K <- NULL
  
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

#' @export
mediation.scan.qr <- function(mediation.qr.object, 
                              data, phenotype,
                              return.allele.effects=FALSE,
                              chr="all", id="SUBJECT.NAME",
                              debug.single.fit=FALSE, use.progress.bar=TRUE,
                              ...){
  model <- mediation.qr.object$model
  subjects <- mediation.qr.object$subjects
  founders <- mediation.qr.object$founders
  num.founders <- length(founders)
  chromatin <- names(mediation.qr.object$qr.alt.list)
  chromatin.chr <- mediation.qr.object$chr
  rh.formula <- mediation.qr.object$formula
  
  if(model == "full" & return.allele.effects){
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available with additive model\n",
        "Setting return.allele.effects to FALSE\n")
  }
  
  rownames(data) <- as.character(data[,id])
  data <- data[subjects,]
  n <- nrow(data)
  
  if(chr[1] != "all"){
    chromatin <- chromatin[chromatin.chr %in% chr]
  }
  
  formula.string <- paste(phenotype, rh.formula)
  formula <- formula(formula.string)
  
  allele.effects <- NULL
  p.vec <- rep(NA, length(chromatin))
  
  if(return.allele.effects){ 
    allele.effects <- matrix(NA, nrow=length(founders), ncol=length(chromatin),
                             dimnames=list(founders, chromatin))
  }
  y <- as.vector(model.frame(formula, data=data)[,1])
  names(y) <- subjects
  # Progress bar
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(chromatin), style=3) }
  for(i in 1:length(chromatin)){
    p.vec[i] <- get.f.stat.p.val(qr.alt=mediation.qr.object$qr.alt.list[[i]], 
                                 qr.null=mediation.qr.object$qr.0.list[[i]], 
                                 y=y)
    if(is.nan(p.vec[i])){ p.vec[i] <- 1 }
    if(return.allele.effects){
      allele.effects[,i] <- get.allele.effects.from.fixef.eQTL(qr.alt=mediation.qr.object$qr.list[[i]], 
                                                               y=y, 
                                                               founders=founders,
                                                               intercept.allele=mediation.qr.object$intercept.allele[i])
    }
    if(debug.single.fit){ browser() }
    # Update progress bar
    if(use.progress.bar){ setTxtProgressBar(pb, i) }
  }
  names(p.vec) <- chromatin
  if(!is.null(mediation.qr.object$condition.loci)){
    formula.string <- paste(Reduce(paste, deparse(formula)), paste(mediation.qr.object$condition.loci, collapse="+"), sep="+")
  }
  output <- list(LOD=NULL,
                 p.value=p.vec,
                 df=NULL,
                 pos=list(Mb=mediation.qr.object$pos$Mb, 
                          cM=mediation.qr.object$pos$cM),
                 loci=chromatin, 
                 chr=chromatin.chr,
                 allele.effects=allele.effects,
                 y=y,
                 formula=formula.string,
                 model.type=model,
                 p.value.method="ANOVA",
                 locus.effect.type="fixed",
                 n=length(y))
  return(output)
}

#' @export
run.qr.permutation.threshold.mediation.scans <- function(perm.ind.matrix, mediation.qr.object, genomecache,
                                                         phenotype, data,
                                                         keep.full.scans=FALSE, scan.index=NULL, id="SUBJECT.NAME",
                                                         chr="all", just.these.loci=NULL, use.progress.bar=TRUE,
                                                         ...){
  
  if(is.null(scan.index)){ scan.index <- 1:ncol(perm.ind.matrix) }
  
  model <- mediation.qr.object$model
  condition.loci <- mediation.qr.object$condition.loci
  locus <- mediation.qr.object$locus
  
  loci <- names(mediation.qr.object$qr.list)
  rh.formula <- mediation.qr.object$formula
  loci.chr <- mediation.qr.object$chr
  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  formula.string <- paste(phenotype, rh.formula, "+", id)
  formula <- formula(formula.string)
  
  full.p <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- matrix(NA, nrow=length(scan.index), ncol=length(loci))
    colnames(full.p) <- loci
    these.pos <- list(Mb=qr.object$pos$Mb[loci],
                      cM=qr.object$pos$cM[loci])
  }
  min.p <- max.p <- rep(NA, length(scan.index))
  
  y <- model.frame(formula, data=data)
  names(y)[1] <- "y"
  y <- y[,c(1, ncol(y))]
  formula.string <- Reduce(paste, deparse(formula))
  perm.formula <- formula(paste0("y ~ ", unlist(strsplit(formula.string, split="~"))[-1]))
  for(i in 1:length(scan.index)){
    ## Permuting all variables but covariates
    permute.var <- !(colnames(data) %in% all.vars(formula)[-1])
    perm.data <- data[perm.ind.matrix[,scan.index[i]],permute.var]
    nonperm.data <- data[,!permute.var]
    perm.data <- cbind(perm.data, nonperm.data)
    
    this.data <- merge(x=y, y=perm.data, by=id, all.x=TRUE)
    this.mediation.qr.object <- extract.mediation.qr(genomecache=genomecache, id=id,
                                                     data=this.data, formula=as.formula(rh.formula), model=model, condition.loci=condition.loci,
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
                                chr=loci.chr, 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=list(min=min.p,
                                               max=max.p))))
}