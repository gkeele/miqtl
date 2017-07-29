#' @export
extract.qr <- function(genomecache, pheno.id="SUBJECT.NAME", geno.id="SUBJECT.NAME",
                       data, formula, model=c("additive", "full"),
                       chr="all", just.these.loci=NULL){
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
  
  pb <- txtProgressBar(min=0, max=length(loci), style=3)
  qr.list <- list()
  for(i in 1:length(loci)){
    X <- h$getLocusMatrix(loci[i], model=model, subjects=subjects)
    keep.col <- 1:ncol(X)
    max.column <- which.max(colSums(X, na.rm=TRUE))[1]
    keep.col <- keep.col[keep.col != max.column]
    X <- cbind(X.0, X[,keep.col])
    qr.list[[i]] <- qr(X)
    setTxtProgressBar(pb, i)
  }
  names(qr.list) <- loci
  
  qr.object <- list(qr.list=qr.list,
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

get.allele.effects.from.fixef.eQTL <- function(qr.alt, y, founders, 
                                               center=TRUE, scale=FALSE){
  regression.effects <- qr.coef(qr=qr.alt, y=y)
  effects <- regression.effects[founders]
  intercept <- regression.effects["(Intercept)"]
  names(effects) <- founders
  
  effects <- effects + intercept
  effects[which(is.na(effects))] <- intercept
  return(as.vector(scale(effects, center=center, scale=scale)))
}

#' @export
scan.qr <- function(qr.object, 
                    data, formula,
                    return.allele.effects=FALSE,
                    chr="all", pheno.id="SUBJECT.NAME", geno.id="SUBJECT.NAME",
                    just.these.loci=NULL,
                    debug.single.fit=FALSE,
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
  
  rownames(data) <- as.character(data[,geno.id])
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
  pb <- txtProgressBar(min=0, max=length(loci), style=3)
  for(i in 1:length(loci)){
    p.vec[i] <- get.f.stat.p.val(qr.alt=qr.object$qr.list[[i]], 
                                 qr.null=qr.object$qr.0, 
                                 y=y, n=n)
    if(return.allele.effects){
      allele.effects[,i] <- get.allele.effects.from.fixef.eQTL(qr.alt=qr.object$qr.list[[i]], y=y, founders=founders)
    }
    if(debug.single.fit){ browser() }
    # Update progress bar
    setTxtProgressBar(pb, i)
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
                 locus.effect.type="fixed")
  return(output)
}



