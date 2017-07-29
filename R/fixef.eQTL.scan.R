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
  data.and.K <- make.processed.data(formula=formula, data=data, 
                                    cache.subjects=cache.subjects, K=K, 
                                    pheno.id=pheno.id, geno.id=geno.id)
  data <- data.and.K$data
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
  
  null.formula <- make.null.formula(formula=formula, do.augment=FALSE)
  X.0 <- model.matrix(null.formula, data=data)
  qr.0 <- qr(X.0)
  
  qr.list <- list()
  for(i in 1:length(loci)){
    X <- h$getLocusMatrix(loci[i], model=model, subjects=cache.subjects)
    keep.col <- 1:ncol(X)
    max.column <- which.max(colSums(X, na.rm=TRUE))[1]
    keep.col <- keep.col[keep.col != max.column]
    X <- X[,keep.col]
    qr.list[[i]] <- qr(X)
  }
  names(qr.list) <- loci
  
  qr.object <- list(qr.list=qr.list,
                    qr.0=qr.0,
                    pos=list(cM=h$getLocusStart(loci, scale="cM"),
                             Mb=h$getLocusStart(loci, scale="Mb")),
                    model=model,
                    founder=h$getFounders())
}

qr.scan <- function(qr.object, 
                    data, formula,
                    return.allele.effects=FALSE,
                    chr="all", pheno.id="SUBJECT.NAME", geno.id="SUBJECT.NAME",
                    just.these.loci=NULL,
                    debug.single.fit=FALSE,
                    ...){
  model <- qr.object$model
  if(model == "full" & return.allele.effects){
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available in additive model", "\n",
        "Setting return.allele.effects to FALSE", "\n")
  }
  
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model="additive"))
  data.and.K <- make.processed.data(formula=formula, data=data, 
                                    cache.subjects=cache.subjects, K=K, 
                                    pheno.id=pheno.id, geno.id=geno.id)
  data <- data.and.K$data
  K <- data.and.K$K
  cache.subjects <- unique(as.character(data[,geno.id]))
  if(!is.null(weights)){ weights <- weights[as.character(data[,pheno.id])] }
  
  loci.chr <- h$getChromOfLocus(loci)
  if(chr[1] != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  augment.indicator <- NULL
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
  original.n <- nrow(data)
  old.data <- data
  
  ###### Augmentation
  if(do.augment){
    augment.n <- ifelse(model=="additive", num.founders, num.founders + choose(num.founders, 2))
    augment.indicator <- c(rep(0, original.n), rep(1, augment.n))
    if(!use.full.null){
      data <- make.simple.augment.data(data=data, formula=formula, augment.n=augment.n)
      data <- data.frame(data, augment.indicator=augment.indicator)
      K <- make.simple.augment.K(K=K, augment.n=augment.n)
    }
    if(use.full.null){
      no.augment.K <- K
      K <- make.full.null.augment.K(K=no.augment.K, original.n=original.n, augment.n=augment.n)
      data <- make.full.null.augment.data(formula=formula, data=data, no.augment.K=no.augment.K, use.par=use.par, brute=brute,
                                          original.n=original.n, augment.n=augment.n, weights=weights)
    }
    weights <- make.augment.weights(data=data, weights=weights, augment.n=augment.n, added.data.points=added.data.points)
  }
  ###### Null model
  ## check for LMER notation
  use.lmer <- check.for.lmer.formula(null.formula)
  
  ####################### CHECKS
  ###### check that full directory has data
  if(model == "full" | use.multi.impute){
    if(!file.exists(paste0(genomecache, "/full/chr", loci.chr[1], "/data/", loci[1], ".RData"))){
      stop("Full model probabilities not available in genome cache, only additive ROP can be fit", call.=FALSE)
    }
  }
  ###### check that full model and return.allele.effects aren't both specified
  if(model == "full" & return.allele.effects){
    return.allele.effects <- FALSE
    cat("Allele effects from regression models currently only available with additive model\n",
        "Setting return.allele.effects to FALSE\n")
  }
  ###### check p-value method
  if(use.lmer & p.value.method == "ANOVA"){
    cat("ANOVA not currently supported with our implementation of LMER, swithcing to LRT\n")
    p.value.method <- "LRT"
  }
  else if(p.value.method == "ANOVA" & (!is.null(K) | !locus.as.fixed)){
    cat("Standard ANOVA F-test not valid with mixed effect model, swithcing to LRT\n")
    p.value.method <- "LRT"
  }
  ###### check that both sparse random effect and kinship random effect are not being specified together
  if(use.lmer & !is.null(K)){
    stop("Cannot use LMER sparse random effects AND a non-sparse random effect", call.=FALSE)
  }
  ###### check that both lmer and random locus effect are not being specified together
  if(use.lmer & !locus.as.fixed){
    stop("Cannot use LMER sparse random effects AND fit locus effect as random", call.=FALSE)
  }
  ###### check that use.fix.par=TRUE if locus.as.fixed=FALSE
  if(!use.fix.par & !locus.as.fixed){
    cat("standard ANOVA F-test not valid with mixed effect model, swithcing to LRT\n")
    use.fix.par <- TRUE
  }
  ###### Switching p-value LRT mode for a random effect
  if(!locus.as.fixed & p.value.method == "LRT"){
    p.value.method <- "LRT.random.locus"
  }
  
  ###### Null model fits
  if(use.lmer){
    fit0 <- lmmbylmer(formula=null.formula, data=data, REML=FALSE, weights=weights)
    fit0.REML <- lmmbylmer(formula=null.formula, data=data, REML=TRUE, weights=weights)
    fix.par <- NULL
  }
  else{
    ## No kinship effect - weights or no weights
    if(is.null(K)){
      fit0 <- lmmbygls(formula=null.formula, data=data, eigen.K=NULL, K=NULL, pheno.id=pheno.id,
                       use.par="h2", fix.par=0, weights=weights, brute=brute)
      fit0.REML <- lmmbygls(formula=null.formula, data=data, eigen.K=NULL, K=NULL, pheno.id=pheno.id,
                            use.par="h2.REML", fix.par=0, weights=weights, brute=brute)
    }
    ## Kinship effect - weights or no weights
    else{
      ###### Handling replicates
      if(pheno.id != geno.id){
        Z <- model.matrix(process.random.formula(geno.id=geno.id), data=data)
        K <- Z %*% K %*% t(Z)
        rownames(K) <- colnames(K) <- as.character(data[,pheno.id])
      }
      ###### Handling constant weights at all loci
      if(!is.null(weights)){
        J <- weights^(1/2) * t(weights^(1/2) * K)
        eigen.J <- process_eigen_decomposition(eigen.decomp=eigen(J))
        fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.J, K=K, use.par=use.par, weights=weights, brute=brute)
        fit0.REML <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.J, K=K, use.par="h2.REML", weights=weights, brute=brute)
      }
      else{
        eigen.K <- process_eigen_decomposition(eigen.decomp=eigen(K))
        fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.K, K=K, use.par=use.par, weights=weights, brute=brute)
        fit0.REML <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=weights, brute=brute)
      }
    }
    ####### EMMA or EMMAX  
    if(use.fix.par){
      fix.par <- ifelse(locus.as.fixed, fit0$h2, fit0.REML$h2)
    }
    if(!use.fix.par){
      fix.par <- NULL
    }
  }
  MI.LOD <- MI.p.value <- allele.effects <- NULL
  LOD.vec <- p.vec <- df <- rep(NA, length(loci))
  null.data <- data
  
  if(return.allele.effects){ 
    if(use.multi.impute){
      allele.effects <- array(NA, dim=c(length(founders), length(loci), num.imp),
                              dimnames=list(founders, loci, paste0("imp", 1:num.imp)))
    }
    else{
      allele.effects <- matrix(NA, nrow=length(founders), ncol=length(loci),
                               dimnames=list(founders, loci))
    }
  }
  ## Prepping link between phenotype and genotype (necessary for imputation in multiple imputations)
  impute.map <- data.frame(data[,pheno.id], data[,geno.id])
  names(impute.map) <- c(pheno.id, geno.id)
  non.augment.subjects <- as.character(data[,geno.id])[grep(pattern="augment", x=as.character(data[,geno.id]), invert=TRUE)]
  
  
  # Progress bar
  pb <- txtProgressBar(min=0, max=length(loci), style=3)
  for(i in 1:length(loci)){
    if(use.multi.impute){
      if(i == 1){ # only at the beginning
        MI.LOD <- MI.p.value <- matrix(NA, nrow=num.imp, ncol=length(loci))
      }
      diplotype.prob.matrix <- h$getLocusMatrix(loci[i], model="full", subjects=non.augment.subjects)
      if(do.augment){
        if(model=="additive"){
          augment.matrix <- matrix(0, nrow=augment.n, ncol=choose(augment.n, 2) + augment.n)
          for(k in 1:augment.n){
            augment.matrix[k, k] <- 1
          }
        }
        if(model=="full"){
          augment.matrix <- diag(augment.n)
        }
        sample.names <- rownames(diplotype.prob.matrix)
        diplotype.prob.matrix <- rbind(diplotype.prob.matrix, augment.matrix)
        rownames(diplotype.prob.matrix) <- c(sample.names, paste0("augment.obs", 1:augment.n))
      }
      if(locus.as.fixed){ fit0.for.mi <- fit0 }
      else{ fit0.for.mi <- fit0.REML }
      fit1 <- multi.imput.lmmbygls(num.imp=num.imp, data=data, formula=formula, weights=weights, locus.as.fixed=locus.as.fixed, return.allele.effects=return.allele.effects,
                                   model=model, p.value.method=p.value.method, founders=founders, diplotype.probs=diplotype.prob.matrix, pheno.id=pheno.id, 
                                   use.lmer=use.lmer, impute.map=impute.map,
                                   use.par=use.par, fix.par=fix.par, fit0=fit0.for.mi, do.augment=do.augment, 
                                   brute=brute, seed=seed) 
      MI.LOD[,i] <- fit1$LOD
      MI.p.value[,i] <- fit1$p.value
      LOD.vec[i] <- median(fit1$LOD)
      p.vec[i] <- median(fit1$p.value)
      
      if(return.allele.effects){
        allele.effects[,i,] <- fit1$allele.effects
      }
    }
    if(!use.multi.impute){
      X <- h$getLocusMatrix(loci[i], model=model, subjects=non.augment.subjects)
      keep.col <- 1:ncol(X)
      if(locus.as.fixed){
        max.column <- which.max(colSums(X, na.rm=TRUE))[1]
        keep.col <- keep.col[keep.col != max.column]
        X <- X[,keep.col]
      }
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
      locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
      if(do.augment){
        X.names <- rownames(X)
        if(model=="additive"){
          X <- rbind(X, 2*diag(augment.n)[,keep.col])
        }
        if(model=="full"){
          X <- rbind(X, diag(augment.n)[,keep.col])
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
      }
      
      data <- cbind(null.data, X)
      if(use.lmer){
        fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
        LOD.vec[i] <- log10(exp(as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))))
        p.vec[i] <- pchisq(q=-2*(as.numeric(logLik(fit0)) - as.numeric(logLik(fit1))), df=length(fixef(fit1))-length(fixef(fit0)), lower.tail=FALSE)
      }
      else{
        if(locus.as.fixed){
          fit1 <- lmmbygls(formula=locus.formula, data=data, pheno.id=pheno.id,
                           eigen.K=fit0$eigen.K, K=fit0$K, weights=weights,
                           use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                           brute=brute)
          LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
          p.vec[i] <- get.p.value(fit0=fit0, fit1=fit1, method=p.value.method)
          df[i] <- fit1$rank
          
          if(return.allele.effects){
            allele.effects[,i] <- get.allele.effects.from.fixef(fit=fit1, founders=founders, allele.in.intercept=founders[max.column])
          }
        }
        else{
          fit1 <- lmmbygls.random(formula=null.formula, data=data, pheno.id=pheno.id,
                                  eigen.K=fit0$eigen.K, K=fit0$K, Z=X, weights=weights,
                                  use.par="h2", null.h2=fix.par,
                                  brute=brute)
          LOD.vec[i] <- log10(exp(fit1$REML.logLik - fit0.REML$REML.logLik))
          p.vec[i] <- get.p.value(fit0=fit0.REML, fit1=fit1, method=p.value.method)
          df[i] <- 1
          if(return.allele.effects){
            allele.effects[,i] <- get.allele.effects.from.ranef(fit=fit1, founders=founders)
          }
        }
      }
    }
    if(debug.single.fit){ browser() }
    # Print out locus fit
    if(print.locus.fit){ cat(paste("locus", i, "out of", length(loci)), "\n") }
    else{
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
  }
  names(LOD.vec) <- names(p.vec) <- names(df) <- loci
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 MI.LOD=MI.LOD,
                 MI.p.value=MI.p.value,
                 df=df,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 fit0.REML=fit0.REML,
                 allele.effects=allele.effects,
                 y=data$y,
                 formula=formula.string,
                 model.type=model,
                 p.value.method=p.value.method,
                 impute.map=impute.map,
                 locus.effect.type=fit1$locus.effect.type)
  if(length(just.these.loci) == 1){ output$fit1 <- fit1 }
  if(pheno.id != geno.id & !is.null(K)){ rownames(Z) <- as.character(data[, pheno.id]); output$Z <- Z }
  return(output)
}



instability.lm.scan <- function(simple.sample.object,
                                genomecache,
                                model=c("additive", "full"),
                                seed=1,
                                use.ROP=TRUE, num.imp=11, chr="all", just.these.loci=NULL, print.locus.fit=TRUE,
                                ...){
  
  y.matrix <- simple.sample.object$y.matrix
  num.scans <- ncol(y.matrix) - 1
  pheno.id <- simple.sample.object$pheno.id
  pheno.data <- simple.sample.object$data
  null.formula <- simple.sample.object$formula
  null.formula.string <- ifelse(is.formula(null.formula), Reduce(paste, deparse(null.formula)), null.formula)
  
  num.samples <- nrow(y.matrix)
  
  h <- DiploprobReader$new(genomecache)
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
  
  full.results <- matrix(NA, nrow=num.scans, ncol=length(loci))
  colnames(full.results) <- loci
  these.chr <- h$getChromOfLocus(loci)
  these.pos <- list(cM=h$getLocusStart(loci, scale="cM"),
                    Mb=h$getLocusStart(loci, scale="Mb"))
  founders <- gsub(pattern="/", replacement=".", x=h$getFounders(), fixed=TRUE)
  ####################### closures (functions defined within function)
  get.qr.null <- function(){
    qr.null <- qr.alt
    num.gen.columns <- ifelse(model=="additive", length(founders)-1, length(founders)+choose(length(founders)-1, 2))
    total.columns <- ncol(design.matrix)
    qr.null$qr <- qr.null$qr[,1:(total.columns-num.gen.columns), drop=FALSE]
    qr.null$rank <- length(1:(total.columns-num.gen.columns))
    qr.null$qraux <- qr.null$qraux[1:(total.columns-num.gen.columns)]
    qr.null$pivot <- qr.null$pivot[1:(total.columns-num.gen.columns)]
    return(qr.null)
  }
  run.locus.fits <- function(){
    get.f.stat.p.val <- function(qr.alt, qr.null, y){
      rss0 <- sum(qr.resid(qr.null, y)^2)
      rss1 <- sum(qr.resid(qr.alt, y)^2)
      df1 <- qr.alt$rank - qr.null$rank
      df2 <- num.samples - qr.alt$rank
      
      mst <- (rss0 - rss1)/df1
      mse <- rss1/df2
      f.stat <- mst/mse
      p.val <- pf(q=f.stat, df1=df1, df2=df2, lower.tail=FALSE)
      return(p.val)
    }
    #######################
    this.locus <- rep(NA, num.scans)
    for(sample in 1:num.scans){
      this.y <- y.matrix[,sample]
      
      this.locus[sample] <- get.f.stat.p.val(qr.alt=qr.alt, qr.null=qr.null, y=this.y)
    }
    return(this.locus)
  }
  
  impute.results <- array(NA, dim=c(num.imp, num.scans, length(loci)))
  for(i in 1:length(loci)){
    this.locus.scan <- rep(NA, num.scans)
    if(use.ROP){
      # ROP
      geno.data <- h$getLocusMatrix(loci[i], model=model, subjects=as.character(pheno.data[,pheno.id]))
      geno.names <- gsub(pattern="/", replacement=".", x=colnames(geno.data), fixed=TRUE)
      set.to.intercept <- which.max(colSums(geno.data))
      combined.data <- data.frame(pheno.data, geno.data[,-set.to.intercept])
      alt.formula.string <- paste("~", unlist(strsplit(paste0(null.formula.string, " + ", paste(geno.names[-set.to.intercept], collapse=" + ")), split="~"))[2])
      design.matrix <- model.matrix(formula(alt.formula.string), data=combined.data)
      
      qr.alt <- qr(design.matrix)
      # get qr.null from qr.alt
      if(i == 1){
        qr.null <- get.qr.null()
      }
      this.locus <- run.locus.fits()
      full.results[,i] <- this.locus
      impute.results <- NULL
    }
    if(!use.ROP){
      set.seed(seed)
      prob.geno.data <- h$getLocusMatrix(loci[i], model="full", subjects=as.character(pheno.data[,pheno.id]))
      prob.geno.data[prob.geno.data < 0] <- 0
      for(imp in 1:num.imp){
        geno.names <- colnames(prob.geno.data)
        geno.data <- t(apply(prob.geno.data, 1, function(x) rmultinom(1, 1, x)))
        colnames(geno.data) <- geno.names
        if(model == "additive"){
          full.to.dosages <- straineff.mapping.matrix()
          geno.data <- geno.data %*% full.to.dosages
          colnames(geno.data) <- founders
          geno.names <- founders
        }
        set.to.intercept <- which.max(colSums(geno.data))
        combined.data <- data.frame(pheno.data, geno.data[,-set.to.intercept])
        alt.formula.string <- paste0(null.formula.string, " + ", paste(geno.names[-set.to.intercept], collapse=" + "))
        design.matrix <- model.matrix(formula(alt.formula.string), data=combined.data)
        qr.alt <- qr(design.matrix)
        if(i == 1){ 
          qr.null <- get.qr.null() 
        }
        this.locus <- run.locus.bs.fits()
        impute.results[imp,,i] <- this.locus
      }
      full.results[,i] <- apply(impute.results[,,i, drop=FALSE], 2, function(x) median(x))
    }
    if(print.locus.fit){
      cat(paste0("Finished locus: ", i, " of ", length(loci), "\n"))
    }
  }
  return(list(impute.results=impute.results, full.results=list(p.values=full.results, chr=these.chr, pos=these.pos),
              formula=null.formula,
              model=model))
}