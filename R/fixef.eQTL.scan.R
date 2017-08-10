#' @export
extract.qr <- function(genomecache, pheno.id="SUBJECT.NAME", geno.id="SUBJECT.NAME",
                       data, formula, model=c("additive", "full"),
                       chr="all", just.these.loci=NULL, use.progress.bar=TRUE){
  K <- NULL
  
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model="additive"))
  cache.subjects <- unique(as.character(data[,geno.id]))
  
  loci.chr <- h$getChromOfLocus(loci)
  if(chr[1] != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  subjects <- as.character(data[,geno.id])
  X.0 <- model.matrix(formula, data=data)
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
                    qr.0=qr.0,
                    chr=loci.chr,
                    pos=list(cM=h$getLocusStart(loci, scale="cM"),
                             Mb=h$getLocusStart(loci, scale="Mb")),
                    model=model,
                    founders=h$getFounders(),
                    subjects=subjects)
}

get.f.stat.p.val <- function(qr.alt, qr.null, y, n){
  rss0 <- sum(qr.resid(qr.null, y)^2)
  rss1 <- sum(qr.resid(qr.alt, y)^2)
  df1 <- qr.alt$rank - qr.null$rank
  df2 <- n - qr.alt$rank
  
  mst <- (rss0 - rss1)/df1
  mse <- rss1/df2
  f.stat <- mst/mse
  p.val <- pf(q=f.stat, df1=df1, df2=df2, lower.tail=FALSE)
  return(p.val)
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

#' @export
scan.qr <- function(qr.object, 
                    data, formula,
                    return.allele.effects=FALSE,
                    chr="all", id="SUBJECT.NAME",
                    just.these.loci=NULL,
                    debug.single.fit=FALSE, use.progress.bar=TRUE,
                    ...){
  model <- qr.object$model
  subjects <- qr.object$subjects
  founders <- qr.object$founders
  num.founders <- length(founders)
  loci <- names(qr.object$qr.list)
  loci.chr <- qr.object$chr

  if(model == "full" & return.allele.effects){
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available with additive model\n",
        "Setting return.allele.effects to FALSE\n")
  }
  
  rownames(data) <- as.character(data[,id])
  data <- data[subjects,]
  n <- nrow(data)

  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  formula.string <- Reduce(paste, deparse(formula))

  allele.effects <- NULL
  p.vec <- rep(NA, length(loci))

  if(return.allele.effects){ 
    allele.effects <- matrix(NA, nrow=length(founders), ncol=length(loci),
                             dimnames=list(founders, loci))
  }
  y <- model.frame(formula, data=data)[,1]
  names(y) <- subjects
  # Progress bar
  if(use.progress.bar){ pb <- txtProgressBar(min=0, max=length(loci), style=3) }
  for(i in 1:length(loci)){
    p.vec[i] <- get.f.stat.p.val(qr.alt=qr.object$qr.list[[i]], 
                                 qr.null=qr.object$qr.0, 
                                 y=y, n=n)
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

#' @export
generate.qr.permutation.index.matrix <- function(qr.scan.object, num.samples, seed=1){
  n <- length(qr.scan.object$y)
  
  set.seed(seed)
  perm.ind.matrix <- replicate(n=num.samples, sample(1:n, replace=FALSE))
  colnames(perm.ind.matrix) <- paste0("perm.", 1:num.samples)
  rownames(perm.ind.matrix) <- names(qr.scan.object$y)
  return(perm.ind.matrix)
}

#' @export
run.qr.permutation.threshold.scans <- function(perm.ind.matrix, qr.object,
                                               keep.full.scans=FALSE, scan.index=NULL, id="SUBJECT.NAME",
                                               formula, data, model=c("additive", "full"),
                                               chr="all", just.these.loci=NULL, use.progress.bar=TRUE,
                                               ...){
  model <- model[1]

  if(is.null(scan.index)){ scan.index <- 1:ncol(perm.ind.matrix) }
  
  loci <- names(qr.object$qr.list)
  loci.chr <- qr.object$chr
  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  full.p <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- matrix(NA, nrow=length(scan.index), ncol=length(loci))
    colnames(full.p) <- loci
    these.pos <- list(Mb=qr.object$pos$Mb[loci],
                      cM=qr.object$pos$cM[loci])
  }
  min.p <- rep(NA, length(scan.index))
  
  formula.string <- Reduce(paste, deparse(formula))
  perm.formula <- formula(paste0("new_y ~ ", unlist(strsplit(formula.string, split="~"))[-1]))
  for(i in scan.index){
    new.y <- data.frame(perm.ind.matrix[,i], rownames(perm.ind.matrix))
    names(new.y) <- c("new_y", id)
    this.data <- merge(x=new.y, y=data, by=id, all.x=TRUE)
    ## Matrix of permutation indexes
    this.data[,all.vars(formula)[1]] <- this.data[,all.vars(formula)[1]][this.data$new_y]
    
    this.scan <- scan.qr(qr.object=qr.object, data=this.data, 
                         formula=perm.formula, model=model,
                         id=id, chr=chr, return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                         ...)
    if(keep.full.scans){
      full.p[i,] <- this.scan$p.value
    }
    min.p[i] <-  min(this.scan$p.value)
    cat("\n", "Threshold scan:", i, "complete", "\n")
  }
  return(list(full.results=list(LOD=NULL,
                                p.value=full.p,
                                chr=loci.chr, 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=min.p)))
}


