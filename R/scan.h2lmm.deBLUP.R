scan.h2lmm.deBLUP <- function(genomecache, data, formula, K,
                              model=c("additive", "full", "diplolasso"),
                              use.par="h2", use.multi.impute=TRUE, num.imp=10, chr="all", brute=TRUE, 
                              seed=1, 
                              weights=NULL, do.augment=FALSE, use.augment.weights=FALSE, use.full.null=FALSE, added.data.points=1, 
                              just.these.loci=NULL,
                              print.locus.fit=FALSE,
                              ...){
  model <- model[1]

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
  null.formula <- make.null.formula()
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
  if(is.null(weights)){
    eigen.K <- eigen(K)
    fit0 <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par=use.par, weights=weights, calc.h2.range.by=calc.h2.range.by, brute=brute, use.chol=FALSE)
    fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=weights, calc.h2.range.by=calc.h2.range.by, brute=brute, use.chol=FALSE)
  }
  if(!is.null(weights)){
    fit0 <- lmmbygls(null.formula, data=data, K=K, use.par=use.par, weights=weights, calc.h2.range.by=calc.h2.range.by, brute=brute, use.chol=TRUE)
    fit0.REML <- NULL
    if(!use.chol){
      fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=weights, calc.h2.range.by=calc.h2.range.by, brute=brute, use.chol=FALSE)
    }
  }
  
  if(is.null(weights)){
    D <- diag(nrow(K)) 
  }
  if(!is.null(weights)){
    D <- diag(weights)
  }
  X <- fit0$x
  Sigma <- K*fit0$tau2.mle + D*fit0$sigma2.mle
  inv.Sigma <- solve(Sigma)
  u <- (K*fit0$tau2.mle) %*% inv.Sigma %*% (diag(nrow(K)) - X %*% solve(t(X) %*% inv.Sigma %*% X) %*% t(X) %*% inv.Sigma) %*% fit0$y
  data$y <- data$y - u
  fix.par <- 0
  M <- diag(weights^-0.5)
  logDetV <- sum(log(weights))
  fit0 <- lmmbygls(null.formula, data=data, K=K, use.par=use.par, fix.par=fix.par, weights=weights, brute=brute, 
                   M=M, logDetV=logDetV)

  LOD.vec <- LOD.median.vec <- LOD.scale.vec <- p.vec <- mi.p.vec <- p.scale.vec <- rep(0, length(loci))
  LOD.mat <- p.mat <- matrix(0, nrow=num.imp, ncol=length(loci))
  null.data <- data
  
  for(i in 1:length(loci)){
    if(use.multi.impute){
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
      fit1 <- multi.imput.lmmbygls(num.imp=num.imp, data=data, formula=formula,
                                   diplotypes=diplotype.matrix, 
                                   model=model, use.par=use.par, fix.par=fix.par, fit0=fit0, do.augment=do.augment, 
                                   brute=brute, seed=seed, weights=weights, use.chol=use.chol, calc.leverage=calc.leverage, use.CC=use.CC) 
      fit1$h2.grid <- NULL
      h2.record[,i] <- fit1$h2
      LOD.vec[i] <- fit1$mean.LOD
      LOD.median.vec[i] <- fit1$median.LOD
      LOD.scale.vec[i] <- fit1$mean.LOD.Lik
      LOD.mat[,i] <- fit1$LOD
      p.vec[i] <- median(fit1$p.value)
      p.scale.vec[i] <- fit1$mean.p.value.scale
      p.mat[,i] <- fit1$p.value
      mi.p.vec[i] <- fit1$imp.mi.p.value
      if(calc.leverage){
        leverage.mat[[i]] <- fit1$leverage
        names(leverage.mat)[i] <- loci[i]
        MI.X[[i]] <- fit1$imp.X
        names(MI.X)[i] <- loci[i]
      }
    }
    if(!use.multi.impute){
      if(model %in% c("additive", "full")){
        X.check <- h$getLocusMatrix(loci[i], model=model)
        X <- h$getLocusMatrix(loci[i], model=model, subjects=data$SUBJECT.NAME[1:original.n])
        max.column <- which.max(colSums(X, na.rm=TRUE))[1]
        X <- X[,-max.column]
      }
      if(model == "diplolasso"){
        X.dosages <- h$getLocusMatrix(loci[i], model="additive", subjects=data$SUBJECT.NAME[1:original.n])
        max.column <- which.max(colSums(X.dosages))[1]
        X <- cbind(X.dosages[,-max.column],
                   h$getLocusMatrix(loci[i], model="full", subjects=data$SUBJECT.NAME[1:original.n])[,-(1:num.founders)])
      }
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
      locus.formula <- make.alt.formula()
      if(do.augment){
        X.names <- rownames(X)
        if(model=="additive"){
          X <- rbind(X, 2*diag(augment.n)[,-max.column])
        }
        if(model=="full"){
          X <- rbind(X, diag(augment.n)[,-max.column])
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
      }
      if(!do.augment){
        if(use.SIC.weights){
          weights <- 1/loci.SIC[,i]
        }
      }
      data <- cbind(null.data, X)
      use.diplolasso <- FALSE
      diplolasso.penalty.factor <- NULL
      if (model == "diplolasso"){
        diplolasso.penalty.factor <- c(rep(0, ncol(fit0$x) + num.founders - 1), rep(1, choose(num.founders, 2)))
        use.diplolasso <- TRUE
      }
      ## final check for any NAs
      data <- data[!apply(data, 1, function(x) any(is.na(x))),]
      K <- K[as.character(data$SUBJECT.NAME), as.character(data$SUBJECT.NAME)]
      fit1 <- lmmbygls(locus.formula, data=data, 
                       eigen.K=fit0$eigen.K, K=fit0$K, 
                       use.diplolasso=use.diplolasso, penalty.factor=diplolasso.penalty.factor,
                       use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                       brute=brute, 
                       weights=weights, calc.h2.range.by=calc.h2.range.by, use.chol=use.chol)
      LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
      p.vec[i] <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
    }
    if(print.locus.fit){
      cat("locus", i, "\n")
    }
  }
  output <- list(LOD=LOD.vec,
                 LOD.median=LOD.median.vec,
                 LOD.scale=LOD.scale.vec,
                 mi.LOD=LOD.mat,
                 p.value=p.vec,
                 p.value.scale=p.scale.vec,
                 mi.p.mat=p.mat,
                 mi.p.vec=mi.p.vec,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 fit0.REML=fit0.REML,
                 leverage=leverage.mat,
                 y=fit1$y,
                 formula=formula.string,
                 model.type=model)
  return(output)
}

