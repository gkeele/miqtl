simulate.from.average.over.imputations <- function(){
  for(i in 1:num.av.over){
    if(model=="additive"){
      X <- t(apply(P, 1, function(x) rmultinom(1, 1, x))) %*% full.to.dosages
      if(i == 1){ max.column <- which.max(colSums(X))[1] }
      X <- X[,-max.column]
      colnames(X) <- h$getFounders()[-max.column]
    }
    if(model=="full"){
      X <- t(apply(P, 1, function(x) rmultinom(1, 1, x)))
      if(i == 1){ max.column <- which.max(colSums(X))[1] }
      X <- X[,-max.column]
      colnames(X) <- colnames(P)[-max.column]
    }
    make.alt.formula <- function(){
      this.formula.string <- Reduce(paste, deparse(formula))
      this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
      this.formula.string <- ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string)
      this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + "), sep=" + "))
      return(this.formula)
    }
    locus.formula <- make.alt.formula()
    use.data <- cbind(data, X)
    fit1 <- lmmbygls(locus.formula, data=use.data, eigen.K=eigen.K, K=K, 
                     use.par=use.par, fix.par=fix.par, use.diplolasso=FALSE,
                     brute=FALSE, weights=NULL, use.chol=FALSE)
    cat("\t")
    cat(i, " ")
    keep <- !is.na(fit1$coefficients)
    q.matrix[,i] <- fit1$x[,keep] %*% fit1$coefficients[keep]
    browser()
    sigma2.vec[i] <- fit1$sigma2.mle
    tau2.vec[i] <- fit1$tau2.mle
  }
  cat("\n")
  q.bar <- rowMeans(q.matrix)
  sigma2.bar <- mean(sigma2.vec)
  tau2.bar <- mean(tau2.vec)
  
  y.bs.matrix <- matrix(NA, nrow=nrow(data), ncol=num.bs.samples)
  for(i in 1:num.bs.samples){
    u <- t(chol.K) %*% rnorm(n=nrow(data), mean=0, sd=sqrt(tau2.bar))
    e <- rnorm(n=nrow(data), mean=0, sd=sqrt(sigma2.bar))
    y.bs.matrix[,i] <- q.bar + u + e
  }
  colnames(y.bs.matrix) <- paste0(bs.phenotype.prefix, 1:num.bs.samples)
  return(y.bs.matrix)
}

regress.out.average.qtl <- function(){
  set.seed(seed)
  qtl.pred.matrix <- matrix(NA, nrow=nrow(this.data), ncol=num.av.over)
  for(i in 1:num.av.over){
    if(model=="additive"){
      X <- t(apply(P, 1, function(x) rmultinom(1, 1, x))) %*% full.to.dosages
      max.column <- which.max(colSums(X))[1]
      X <- X[,-max.column]
      colnames(X) <- h$getFounders()[-max.column]
    }
    if(model=="full"){
      X <- t(apply(P, 1, function(x) rmultinom(1, 1, x)))
      max.column <- which.max(colSums(X))[1]
      X <- X[,-max.column]
      colnames(X) <- colnames(P)[-max.column]
    }
    make.alt.formula <- function(){
      this.formula.string <- Reduce(paste, deparse(formula))
      this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
      this.formula.string <- ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string)
      this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + "), sep=" + "))
      return(this.formula)
    }
    locus.formula <- make.alt.formula()
    use.data <- cbind(this.data, X)
    #browser()
    fit1 <- lmmbygls(locus.formula, data=use.data, eigen.K=NULL, K=this.K, 
                     use.par=use.par, fix.par=fix.par, use.diplolasso=FALSE,
                     brute=FALSE, weights=NULL, use.chol=FALSE)
    cat("\t")
    cat(i, " ")
    keep <- names(fit1$coefficients[colnames(X)])[!is.na(fit1$coefficients[colnames(X)])]
    qtl.pred.matrix[,i] <- fit1$y - fit1$x[,keep] %*% fit1$coefficients[keep]
  }
  cat("\n")
  qtl.pred.bar <- rowMeans(qtl.pred.matrix)
  return(qtl.pred.bar)
}