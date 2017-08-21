gls.fit <- function(X, y, M, logDetV, ...){
  n  <- nrow(X)
  # Tranform to uncorrelated Sigma form
  if(!is.null(M)){
    MX <- M %*% X
    My <- M %*% y
  }
  if(is.null(M)){ # more efficient than MX <- diag(nrow(X)) %*% X
    MX <- X
    My <- y
  }
  
  fit <- lm.fit(MX, My, ...)
  fit$rss <- sum(fit$residuals^2)
  fit$sigma2 <- fit$rss/fit$df.residual
  fit$sigma2.mle <- fit$rss/n
  fit$uncorr.logLik <- -0.5*n*log(2*pi) - 0.5*n*log(fit$sigma2.mle) - 0.5*n
  # Correct logLik to value under correlated Sigma
  fit$logLik <- fit$uncorr.logLik - 0.5*logDetV
  fit$logDetV <- logDetV
  fit$M <- M
  return(fit)
}

#' @export
lmmbygls <- function(formula, data=NULL, 
                     y=NULL, X=NULL,
                     K=NULL, eigen.K=NULL, fix.par=NULL,
                     M=NULL, logDetV=NULL, weights=NULL, pheno.id="SUBJECT.NAME",
                     use.par=c("h2", "h2.REML"), 
                     brute=TRUE,
                     subset, na.action,
                     method = "qr",
                     model = TRUE, 
                     contrasts = NULL,
                     verbose = FALSE,
                     ...) 
{
  if(!is.null(data) & is.null(y)){
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$y <- m$X <- m$K <- m$eigen.K <- m$fix.par <- m$use.par <- NULL
    m$M <- m$logDetV <- m$weights <- m$pheno.id <- NULL
    m$method <- m$model <- m$contrasts <- m$verbose <- NULL
    m$brute <- m$... <- NULL
    
    m[[1L]] <- quote(stats::model.frame)
    m <- eval.parent(m)
    if(method == "model.frame"){
      return(m)    
    }
    Terms <- attr(m, "terms")
    
    y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    ids <- data[,pheno.id]
  }
  else{
    ids <- rownames(X)
  }
  n <- nrow(X)
  q <- ncol(X)

  if(is.null(fix.par)){
    if(is.null(K)){ ## No kinship effect setting: K - NULL, eigen.K - NULL
      fix.par <- 0
    }
    else{
      if(is.null(eigen.K)){
        eigen.K <- eigen(K)
      }
      Ut <- t(eigen.K$vectors) # a small optimization
    }
  }
  ## Define local objective function/closure for Brent's optimization
  ### Optimize functions
  h2.fit <- function(h2, logLik.only=TRUE, verbose=FALSE, ...){
    if(is.null(fix.par)){ # EMMA or first time optimization
      d <- h2*eigen.K$values + (1-h2)
      if(is.null(weights)){
        M <- d^-0.5 * Ut 
        logDetV <- sum(log(d))
      }
      else{
        M <- d^-0.5 * t(sqrt(weights) * t(Ut))
        logDetV <- 2*sum(log(1/sqrt(weights))) + sum(log(d)) # maybe right
      }
    }
    else{
      if(fix.par == 0 & is.null(weights)){
        M <- NULL # more efficient than I
        logDetV <- 0
      }
      if(fix.par == 0 & !is.null(weights)){
        M <- diag(sqrt(weights))
        logDetV <- sum(log(1/weights))
      }
    }
    fit <- gls.fit(X=X, y=y, M=M, logDetV=logDetV, ...)
    if(logLik.only){
      if(verbose){
        cat(sep="", "h2 = ", h2, " : logLik = ", fit$logLik, "\n")
      }
      return(fit$logLik)
    }
    fit$h2 <- h2
    fit$M <- M
    fit$logDetV <- logDetV
    return(fit)
  }
  h2.fit.REML <- function(h2, logLik.only=TRUE, verbose=FALSE, ...){
    if(is.null(fix.par)){ # EMMA or first time optimization
      d <- h2*eigen.K$values + (1-h2)
      if(is.null(weights)){
        M <- d^-0.5 * Ut 
        logDetV <- sum(log(d))
      }
      else{
        M <- d^-0.5 * t(sqrt(weights) * t(Ut))
        logDetV <- 2*sum(log(1/sqrt(weights))) + sum(log(d)) # maybe right
      }
    }
    else{
      if(fix.par == 0 & is.null(weights)){
        M <- NULL # more efficient than I
        logDetV <- 0
        d <- rep(1, n)
        Ut <- diag(d)
      }
      if(fix.par == 0 & !is.null(weights)){
        M <- diag(sqrt(weights))
        logDetV <- sum(log(1/weights))
        d <- weights
        Ut <- diag(length(d))
      }
    }
    fit <- gls.fit(X=X, y=y, M=M, logDetV=logDetV, ...)
    
    df <- fit$df.residual
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.reml <- fit$rss/df
    
    ## Check to make sure design matrix is full rank - for REML estimation
    col.keep <- !is.na(fit$coefficients)
    X <- X[,col.keep]
    
    logDetXtX <- log(det(t(X) %*% X))
    if(is.null(weights)){
      logDetXtVinvX <- log(det(t(X)%*%t(Ut)%*%((1/d)*Ut)%*%X))
      fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*logDetXtX - 0.5*logDetXtVinvX - 0.5*logDetV
    }
    else{
      logDetXtVinvX <- log(det(t(X)%*%(sqrt(weights)*t(Ut))%*%((1/d)*Ut)%*%(sqrt(weights)*X)))
      fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*logDetXtX - 0.5*logDetXtVinvX - 0.5*logDetV
    }

    if(logLik.only){
      if(verbose){
        cat(sep="", "h2 = ", h2, " : logLik = ", fit$REML.logLik, "\n")
      }
      return(fit$REML.logLik)
    }
    fit$h2 <- h2
    return(fit)
  }
  ## Optimize using objective function, or otherwise given value of h2
  fit <- NULL
  if(is.null(fix.par)){
    if(verbose){
      cat("Optimizing logLik for parameter\n")
    }
    if(use.par[1] == "h2"){
      peak <- optimize(f=h2.fit, logLik.only=TRUE, verbose=verbose, ..., interval=c(0,1), maximum=TRUE)
      if(brute){
        fit.h2.0 <- h2.fit(h2=0, logLik.only=FALSE, verbose=FALSE)
        if(peak$objective < fit.h2.0$logLik){
          fit <- fit.h2.0
        }
        else{
          fit  <- h2.fit(h2=peak$maximum, logLik.only=FALSE, verbose=FALSE)
        }
      }
      else{
        fit  <- h2.fit(h2=peak$maximum, logLik.only=FALSE, verbose=FALSE)
      }
    }
    if(use.par[1] == "h2.REML"){
      peak <- optimize(f=h2.fit.REML, logLik.only=TRUE, verbose=verbose, ..., interval=c(0,1), maximum=TRUE)
      if(brute){
        fit.h2.0 <- h2.fit.REML(h2=0, logLik.only=FALSE, verbose=FALSE)
        if(peak$objective < fit.h2.0$logLik){
          fit <- fit.h2.0
        }
        else{
          fit <- h2.fit.REML(h2=peak$maximum, logLik.only=FALSE, verbose=FALSE)
        }
      }
      else{
        fit <- h2.fit.REML(h2=peak$maximum, logLik.only=FALSE, verbose=FALSE)
      }
    }
    fit$h2.optimized <- TRUE
  } 
  else{
    if(use.par[1] == "h2"){
      fit <- h2.fit(h2=fix.par, logLik.only=FALSE, verbose=FALSE)
    }
    else if(use.par[1] == "h2.REML"){
      fit <- h2.fit.REML(h2=fix.par, logLik.only=FALSE, verbose=FALSE)
    }
    fit$h2.optimized <- TRUE
  }
  fit$gls.sigma2.mle <- fit$sigma2.mle
  fit$gls.sigma2     <- fit$sigma2
  fit$sigma2.mle <- (1 - fit$h2)*fit$gls.sigma2.mle
  fit$tau2.mle <- fit$h2*fit$gls.sigma2.mle

  if(!is.null(data)){ ## Not within a scan
    fit$terms <- Terms
    fit$call <- call
    if(model){
      fit$model <- m
    }
    fit$na.action <- attr(m, "na.action")
    fit$xlevels <- .getXlevels(Terms, m)
    fit$contrasts <- attr(X, "contrasts")
  }

  fit$locus.effect.type <- "fixed"
  names(y) <- rownames(X) <- ids
  fit$weights <- weights
  fit$x <- X
  fit$y <- y
  fit$eigen.K <- eigen.K
  fit$K <- K
  class(fit) <- "lmmbygls"
  return(fit)
}

## Wrapper for lmmbygls when you have replicate and just little K (genome by genome)
## Wrote this for Yanwei to use with his merge analysis
#' @export
lmmbygls.replicates <- function(formula, data, little.K, pheno.id="SUBJECT.NAME", geno.id, ...){
  use.data <- data[data[,geno.id] %in% colnames(K),]
  use.K <- K[colnames(K) %in% data[,geno.id], colnames(K) %in% data[,geno.id]]
  
  # Making n x n K
  use.data[,geno.id] <- as.character(use.data[,geno.id])
  Z <- model.matrix(process.random.formula(geno.id=geno.id), data=use.data)
  big.K <- Z %*% use.K %*% t(Z)
  rownames(big.K) <- colnames(big.K) <- as.character(use.data[,pheno.id])
  lmmbygls.fit <- lmmbygls(formula=formula, data=use.data, K=big.K,
                           pheno.id=pheno.id, ...)
  return(lmmbygls.fit)
}



