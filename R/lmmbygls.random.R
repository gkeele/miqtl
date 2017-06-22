lmmbygls.random <- function(formula, data, K=NULL, eigen.K=NULL, Z, null.h2,
                            weights=NULL, 
                            use.par="h2",
                            pheno.id="SUBJECT.NAME", brute=TRUE,
                            subset, na.action,
                            method="qr",
                            model=TRUE, 
                            contrasts=NULL,
                            verbose=FALSE,
                            ...) 
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$W <- m$Z <- m$K <- m$eigen.K <- m$null.h2 <- m$use.par <- NULL
  m$method <- m$model <- m$weights <- m$pheno.id <- m$brute <- m$x <- m$y <- m$contrasts <- NULL
  m$... <- NULL
  
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  if(method == "model.frame"){
    return(m)    
  }
  Terms <- attr(m, "terms")
  y <- model.response(m)
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)
  q <- ncol(X)
  
  ids <- data[,pheno.id]
  
  if(is.null(K)){ ## No kinship effect setting: K - NULL, eigen.K - NULL
    d <- rep(1, n)
    Ut <- diag(d)
  }
  else{
    if(is.null(eigen.K)){
      eigen.K <- eigen(K)
    }
    Ut <- t(eigen.K$vectors) # a small optimization
    d <- eigen.K$values
  }
  poly.K <- K
  original.y <- y
  original.X <- X
  original.Z <- Z
  ## Rotations
  if(!is.null(weights)){
    rotate <- t(sqrt(weights) * t(Ut))
    un.rotate <- 1/sqrt(weights)*Ut
  }
  else{
    rotate <- Ut
    #un.rotate <- t(Ut)
  }
  y <- rotate %*% y
  X <- rotate %*% X
  Z <- rotate %*% Z
  K <- Z %*% t(Z)

  ## Define local objective function/closure for Brent's optimization
  ### Optimize functions
  h2.fit <- function(h2, logLik.only=TRUE, verbose=FALSE, ...){
    H <- K*h2 + diag(d*null.h2*(1 - h2) + 1 - null.h2*(1 - h2) - h2)
    chol.H <- chol(H)
    M <- t(solve(chol.H))
    fit <- gls.fit(X=X, y=y, M=M, logDetV=0, ...)
    # REML likelihood
    df <- fit$df.residual
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.reml <- fit$rss/df
    ## Check to make sure design matrix is full rank - for REML estimation
    col.keep <- !is.na(fit$coefficients)
    X <- X[,col.keep]
    #raw.H <- un.rotate %*% H %*% t(un.rotate)
    #chol.raw.H <- chol(raw.H)
    if(is.null(weights)){
      fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*log(det(t(X) %*% X)) - 0.5*log(det(t(X) %*% chol2inv(chol.H) %*% X)) - sum(log(diag(chol.H)))
    }
    else{
      fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*log(det(t(original.X) %*% original.X)) - 0.5*log(det(t(X) %*% un.rotate %*% chol2inv(chol.H) %*% t(un.rotate) %*% X)) - 0.5*log(det(un.rotate %*% H %*% t(un.rotate)))
    }
    if(logLik.only){
      if(verbose){
        cat(sep="", "h2 = ", h2, " : logLik = ", fit$REML.logLik, "\n")
      }
      return(fit$REML.logLik)
    }
    fit$h2 <- h2
    fit$lambda <- h2/(1 - h2)
    return(fit)
  }
  ## Optimize using objective function, or otherwise given value of h2
  fit <- NULL
  if(use.par == "h2"){
    peak <- optimize(f=h2.fit, logLik.only=TRUE, ..., interval=c(0, 1-null.h2), maximum=TRUE)
    browser()
    if(brute){
      fit.h2.0 <- h2.fit(h2=0, logLik.only=FALSE, verbose=FALSE)
      if(peak$objective < fit.h2.0$REML.logLik){
        fit <- fit.h2.0
      }
      else{
        fit <- h2.fit(h2=peak$maximum, logLik.only=FALSE)
      }
    }
    else{
      fit <- h2.fit(h2=peak$maximum, logLik.only=FALSE)
    }
  }
  fit$h2.optimized <- TRUE
  fit$terms <- Terms
  fit$call <- call
  if(model){
    fit$model <- m
  }
  names(y) <- rownames(X) <- rownames(Z) <- ids
  fit$locus.h2 <- fit$h2
  fit$h2 <- null.h2
  fit$locus.effect.type <- "random"
  fit$na.action <- attr(m, "na.action")
  fit$weights <- weights
  fit$x <- original.X
  fit$z <- original.Z
  fit$y <- original.y
  fit$eigen.K <- eigen.K
  fit$K <- poly.K
  fit$xlevels <- .getXlevels(Terms, m)
  fit$contrasts <- attr(X, "contrasts")
  class(fit) <- "lmmbygls"
  return(fit)
}






