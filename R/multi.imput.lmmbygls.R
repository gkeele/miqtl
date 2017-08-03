#' @export
multi.imput.lmmbygls <- function(formula, data, pheno.id="SUBJECT.NAME",
                                 y=NULL, fit0=NULL,
                                 num.imp, founders=founders, X.probs, 
                                 K=NULL, return.allele.effects=FALSE,
                                 use.par, fix.par=NULL, model=c("additive", "full"), p.value.method=c("LRT", "ANOVA"), locus.as.fixed=TRUE,
                                 use.lmer, impute.map,
                                 brute=TRUE, seed=1, do.augment,
                                 weights=NULL){
  model <- model[1]
  p.value.method <- p.value.method[1]
  eigen.K <- logDetV <- M <- allele.effects <- NULL
  if(return.allele.effects){ 
    allele.effects <- matrix(NA, nrow=length(founders), ncol=num.imp,
                             dimnames=list(founders, paste0("imp", 1:num.imp)))
  }
  if(is.null(fit0)){
    null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
    if(use.lmer){
      fit0 <- lmmbylmer(null.formula, data=data, REML=FALSE, weights=weights)
    }
    else{
      fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, K=K, use.par=use.par, brute=brute, weights=weights)
      K <- fit0$K
    }
  }
  if(is.null(weights) & !use.lmer){
    eigen.K <- fit0$eigen.K
  }
  if(!is.null(fix.par) & !use.lmer){
    M <- fit0$M
    logDetV <- fit0$logDetV
  }
  full.to.dosages <- straineff.mapping.matrix()
  
  imp.logLik <- imp.h2 <- imp.df <- imp.LOD <- imp.p.value <- rep(0, num.imp)

  set.seed(seed)
  for(i in 1:num.imp){
    if(model == "additive"){
      if(any(X.probs < 0)){
        X.probs[X.probs < 0] <- 0
        X.probs <- t(apply(X.probs, 1, function(x) x/sum(x)))
      }
      X <- run.imputation(diplotype.probs=X.probs, impute.map=impute.map) %*% full.to.dosages
      colnames(X) <- gsub(pattern="/", replacement=".", x=founders, fixed=TRUE)
    }
    if(model == "full"){
      X <- run.imputation(diplotype.probs=X.probs, impute.map=impute.map)
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X.probs), fixed=TRUE)
    }
    keep.col <- 1:ncol(X)
    if(locus.as.fixed){
      max.column <- which.max(colSums(X, na.rm=TRUE))[1]
      keep.col <- keep.col[keep.col != max.column]
      X <- X[,keep.col]
    }

    locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
    if(use.lmer){
      data <- cbind(null.data, X)
      fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
      imp.logLik[i] <- as.numeric(logLik(fit1))
      imp.h2[i] <- NA
      imp.df[i] <- length(fixef(fit1))
      imp.LOD[i] <- log10(exp(imp.logLik[i] - as.numeric(logLik(fit0))))
      imp.p.value[i] <- pchisq(q=-2*(as.numeric(logLik(fit0)) - imp.logLik), df=imp.df - length(fixef(fit0)), lower.tail=FALSE)
      fit1$locus.effect.type <- "fixed"
    }
    else{
      if(locus.as.fixed){
        X <- cbind(fit0$x, X)
        fit1 <- lmmbygls(formula=locus.formula, pheno.id=pheno.id, eigen.K=eigen.K, K=K,
                         y=y, X=X,
                         logDetV=logDetV, M=M, 
                         use.par="h2", fix.par=fix.par,
                         brute=brute, weights=weights)
        imp.logLik[i] <- fit1$logLik
        imp.h2[i] <- fit1$h2
        imp.df[i] <- fit1$rank
        imp.LOD[i] <- log10(exp(imp.logLik[i] - fit0$logLik))
        imp.p.value[i] <- get.p.value(fit0=fit0, fit1=fit1, method=p.value.method)
        fit1$locus.effect.type <- "fixed"
        
        if(return.allele.effects){
          allele.effects[,i] <- get.allele.effects.from.fixef(fit=fit1, founders=founders, allele.in.intercept=founders[max.column])
        }
      }
      else{
        null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
        fit1 <- lmmbygls.random(formula=null.formula, pheno.id=pheno.id,  
                                y=y, X=fit0$x, K=K, eigen.K=eigen.K, Z=X,
                                use.par="h2", null.h2=fix.par,
                                brute=brute, weights=weights)
        imp.logLik[i] <- fit1$REML.logLik
        imp.h2[i] <- fit1$h2
        imp.df[i] <- fit1$rank
        imp.LOD[i] <- log10(exp(imp.logLik[i] - fit0$REML.logLik))
        imp.p.value[i] <- get.p.value(fit0=fit0, fit1=fit1, method=p.value.method)
        fit1$locus.effect.type <- "random"
        
        if(return.allele.effects){
          allele.effects[,i] <- get.allele.effects.from.ranef(fit=fit1, founders=founders)
        }
      }
    }
  }
  return(list(h2=imp.h2, 
              LOD=imp.LOD,
              p.value=imp.p.value,
              allele.effects=allele.effects,
              locus.effect.type=fit1$locus.effect.type))
}
