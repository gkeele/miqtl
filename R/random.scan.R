#' @export
random.scan.h2lmm <- function(genomecache, data, formula, K,
                              model=c("additive", "full"), use.par=c("h2", "lambda"),
                              use.multi.impute=TRUE, num.imp=10, chr="all", brute=TRUE, 
                              seed=1, 
                              do.augment=FALSE, use.augment.weights=FALSE, use.full.null=FALSE, added.data.points=1, 
                              just.these.loci=NULL,
                              print.locus.fit=FALSE, 
                              ...){
  model <- model[1]
  use.par <- use.par[1]
  
  h <- DiploprobReader$new(genomecache)
  num.founders <- length(h$getFounders())
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data(formula=formula, data=data, cache.subjects=cache.subjects, K=K)
  data <- data.and.K$data
  K <- data.and.K$K
  
  loci <- h$getLoci()
  loci.chr <- h$getChromOfLocus(loci)
  
  if(chr != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  augment.indicator <- NULL
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula, do.augment=do.augment)
  original.n <- nrow(data)
  if(do.augment){
    augment.n <- ifelse(model=="additive", num.founders, num.founders + choose(num.founders, 2))
    augment.indicator <- c(rep(0, original.n), rep(1, augment.n))
    if(!use.full.null){
      data <- make.simple.augment.data(data=data, augment.n=augment.n)
      K <- make.simple.augment.K(K=K, augment.n=augment.n)
    }
    if(use.full.null){
      no.augment.K <- K
      K <- make.full.null.augment.K(K=no.augment.K, original.n=original.n, augment.n=augment.n)
      data <- make.full.null.augment.data(formula=formula, data=data, no.augment.K=no.augment.K, use.par=use.par, brute=brute,
                                          original.n=original.n, augment.n=augment.n, weights=weights)
    }
    if(use.augment.weights){
      weights <- make.augment.weights(data=data, augment.n=augment.n, added.data.points=added.data.points)
    }
  }
  
  eigen.K <- NULL
  if(!is.null(K)){
    eigen.K <- eigen(K)
    fit0 <- lmmbygls.null(null.formula, data=data, K=K, eigen.K=eigen.K, use.par=use.par)
  }
  else{
    fit0 <- lmmbygls.null(null.formula, data=data, K=K, eigen.K=eigen.K, use.par=use.par, fix.par=0)
  }
  fix.null.h2 <- fit0$h2
  fix.null.lambda <- fix.null.h2/(1 - fix.null.h2)

  LOD.median.vec <- LOD.scale.vec <- p.vec <- mi.p.vec <- p.scale.vec <- rep(NA, length(loci))
  LOD.mat <- p.mat <- matrix(0, nrow=num.imp, ncol=length(loci))
  null.data <- data
  
  for(i in 1:length(loci)){
    if(use.multi.impute){
      imp.h2.mat <- matrix(NA, nrow=num.imp, ncol=length(loci))
      diplotype.matrix <- h$getLocusMatrix(loci[i], model="full", subjects=data$SUBJECT.NAME)
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
        sample.names <- rownames(diplotype.matrix)
        diplotype.matrix <- rbind(diplotype.matrix, augment.matrix)
        rownames(diplotype.matrix) <- c(sample.names, paste0("augment.obs", 1:augment.n))
      }
      fit1 <- multi.imput.lmmbygls.random(num.imp=num.imp, data=data, this.formula=null.formula,
                                          diplotypes=diplotype.matrix, 
                                          model=model, use.par=use.par, fit0=fit0, do.augment=do.augment, 
                                          seed=seed) 
      
      h2.record[,i] <- fit1$h2
      p.vec[i] <- median(fit1$p.value)
      p.mat[,i] <- fit1$p.value
    }
    if(!use.multi.impute){
      X <- h$getLocusMatrix(loci[i], model=model, subjects=data$SUBJECT.NAME[1:original.n])
      if(do.augment){
        X.names <- rownames(X)
        if(model=="additive"){
          X <- rbind(X, 2*diag(augment.n))
        }
        if(model=="full"){
          X <- rbind(X, diag(augment.n))
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
      }
      fit1 <- random.lmmbygls(null.formula, data=data,
                              Z=X, fit0=fit0, use.par=use.par,
                              null.h2=fix.null.h2)
      chi.sq <- -2*(fit0$REML.logLik - fit1$REML.logLik)
      p.vec[i] <- ifelse(chi.sq == 0, 1, 0.5*pchisq(q=chi.sq, df=1, lower.tail=FALSE))
    }
    if(print.locus.fit){
      cat("locus", i, "\n")
    }
  }
  output <- list(p.value=p.vec,
                 mi.p.mat=p.mat,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 y=fit1$y,
                 formula=formula.string,
                 model.type=model)
  return(output)
}
