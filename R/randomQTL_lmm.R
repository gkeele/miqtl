
simple.gls.fit <- function(X, y, ...,  M, H, H.inv, logDetH)
{
  n  <- nrow(X)
  # Tranform to uncorrelated Sigma form
  MX <- M %*% X
  My <- M %*% y
  fit <- lm.fit(MX, My, ...)
  return(fit)
}

lmmbygls.null <- function (formula, data, K=NULL, eigen.K=NULL, fix.par=NULL,
                           use.par=c("lambda", "h2"), 
                           subset, na.action,
                           method = "qr",
                           model = TRUE, 
                           contrasts = NULL,
                           ...) 
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$W <- m$K <- m$eigen.K <- m$use.par <- m$fix.par <- NULL
  m$method <- m$model <- m$x <- m$y <- m$contrasts <- NULL
  m$... <- NULL
  
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  if (method == "model.frame"){
    return(m)    
  }
  Terms <- attr(m, "terms")
  y <- model.response(m)
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)

  if(!is.null(K)){
    Ut <- t(eigen.K$vectors)
  }
  ## Define local objective function/closure for Brent's optimization
  ### Optimize functions
  lambda.fit <- function(lambda, logLik.only=TRUE, ...){
    d <- eigen.K$values + lambda
    M <- d^-0.5 * Ut
    logDetH <- sum(log(d))
    H <- t(Ut) %*% (d * Ut)
    H.inv <- t(Ut) %*% ((1/d) * Ut)
    
    fit <- simple.gls.fit(X=X, y=y, M=M, ...)
    # REML likelihood
    df <- fit$df.residual
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.reml <- fit$rss/df
    fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*log(det(t(X) %*% X)) - 0.5*log(det(t(X) %*% H.inv %*% X)) - 0.5*logDetH
    
    if (logLik.only){
      return (fit$REML.logLik)
    }
    fit$lambda <- lambda
    fit$h2 <- lambda/(1 + lambda)
    return(fit)
  }
  h2.fit <- function(h2, logLik.only=TRUE, ...){
    if(!is.null(K)){
      d <- h2*eigen.K$values + (1-h2)
      M <- d^-0.5 * Ut 
    }
    else{
      Ut <- diag(nrow(X))
      d <- rep(1, nrow(X))
      M <- diag(nrow(X))
    }
    fit <- simple.gls.fit(X=X, y=y, M=M, ...)
    # REML likelihood
    df <- fit$df.residual
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.reml <- fit$rss/df
    fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*log(det(t(X) %*% X)) - 0.5*log(det(t(X)%*%t(Ut)%*%diag(1/d)%*%Ut%*%X)) - 0.5*sum(log(d))
    
    if(logLik.only){
      return (fit$REML.logLik)
    }
    fit$h2 <- h2
    fit$lambda <- h2/(1 - h2)
    return(fit)
  }
  
  ## Optimize using objective function, or otherwise given value of h2
  fit <- NULL
  if (use.par[1] == "lambda"){
    if(is.null(fix.par)){
      peak <- optimize(f=lambda.fit, logLik.only=TRUE, ..., interval=c(0, exp(10)), maximum=TRUE)
      fit  <- lambda.fit(lambda=peak$maximum, logLik.only=FALSE)
      fit$h2.optimized <- TRUE
    }
    else{
      fit  <- lambda.fit(lambda=fix.par, logLik.only=FALSE)
    }
  }
  if (use.par[1] == "h2"){
    if(is.null(fix.par)){
      peak <- optimize(f=h2.fit, logLik.only=TRUE, ..., interval=c(0,1), maximum=TRUE)
      fit  <- h2.fit(h2=peak$maximum, logLik.only=FALSE)
      fit$h2.optimized <- TRUE
    }
    else{
      fit  <- h2.fit(h2=fix.par, logLik.only=FALSE)
    }
  }

  fit$terms <- Terms
  fit$call <- call
  if (model){
    fit$model <- m
  }
  fit$na.action <- attr(m, "na.action")
  fit$x <- X
  fit$y <- y
  fit$eigen.K <- eigen.K
  fit$K <- K
  fit$xlevels <- .getXlevels(Terms, m)
  fit$contrasts <- attr(X, "contrasts")
  class(fit) <- "lmmbygls"
  return(fit)
}

random.lmmbygls <- function (formula, data, Z, 
                             null.h2, fit0,
                             use.par,
                             subset, na.action,
                             method = "qr",
                             model = TRUE, 
                             contrasts = NULL,
                             ...) 
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$W <- m$Z <- m$null.h2 <- m$use.par <- m$fit0 <- NULL
  m$method <- m$model <- m$x <- m$y <- m$contrasts <- NULL
  m$... <- NULL
  
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  if (method == "model.frame"){
    return(m)    
  }
  Terms <- attr(m, "terms")
  y <- model.response(m)
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)
  q <- ncol(X)
  
  eigen.K <- fit0$eigen.K
  if(!is.null(eigen.K)){
    Ut <- t(eigen.K$vectors) # a small optimization
    d <- eigen.K$values
  }
  else{
    Ut <- diag(nrow(X))
    d <- rep(1, nrow(X))
  }
  y <- Ut %*% y
  X <- Ut %*% X
  Z <- Ut %*% Z
  K <- Z %*% t(Z)
  
  ## Define local objective function/closure for Brent's optimization
  ### Optimize functions
  lambda.fit <- function(lambda, logLik.only=TRUE, ...){
    null.lambda <- null.h2/(1 - null.h2)
    H <- K*lambda + diag(d*null.lambda + 1)
    chol.H <- chol(H)
    M <- t(solve(chol.H))
    fit <- simple.gls.fit(X=X, y=y, M=M, ...)
    # REML likelihood
    df <- fit$df.residual
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.reml <- fit$rss/df
    fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*log(det(t(X) %*% X)) - 0.5*log(det(t(X) %*% chol2inv(chol.H) %*% X)) - 0.5*2*sum(log(diag(chol.H)))
    if (logLik.only){
      return (fit$REML.logLik)
    }
    fit$lambda <- lambda
    fit$h2 <- lambda/(1 + lambda)
    return(fit)
  }
  h2.fit <- function(h2, logLik.only=TRUE, ...){
    H <- K*h2 + diag(d*null.h2*(1 - h2) + 1 - null.h2*(1 - h2) - h2)
    chol.H <- chol(H)
    M <- t(solve(chol.H))
    fit <- simple.gls.fit(X=X, y=y, M=M, ...)
    # REML likelihood
    df <- fit$df.residual
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.reml <- fit$rss/df
    fit$REML.logLik <- -(0.5*df)*(log(2*pi) + log(fit$sigma2.reml) + 1) + 0.5*log(det(t(X) %*% X)) - 0.5*log(det(t(X) %*% chol2inv(chol.H) %*% X)) - 0.5*2*sum(log(diag(chol.H)))
    if (logLik.only){
      return (fit$REML.logLik)
    }
    fit$h2 <- h2
    fit$lambda <- h2/(1 - h2)
    return(fit)
  }
  ## Optimize using objective function, or otherwise given value of h2
  fit <- NULL
  if(use.par == "lambda"){
    peak <- optimize(f=lambda.fit, logLik.only=TRUE, ..., interval=c(0, exp(10)), maximum=TRUE)
    if(fit0$REML.logLik > peak$objective){
      fit <- fit0
    }
    else{
      fit <- lambda.fit(lambda=peak$maximum, logLik.only=FALSE)
    }
  }
  if(use.par == "h2"){
    peak <- optimize(f=h2.fit, logLik.only=TRUE, ..., interval=c(0, 1), maximum=TRUE)
    if(fit0$REML.logLik > peak$objective){
      fit <- fit0
    }
    else{
      fit <- h2.fit(h2=peak$maximum, logLik.only=FALSE)
    }
  }
  fit$h2.optimized <- TRUE
  fit$terms <- Terms
  fit$call <- call
  if (model){
    fit$model <- m
  }
  fit$na.action <- attr(m, "na.action")
  fit$x <- X
  fit$y <- y
  fit$eigen.K <- eigen.K
  fit$xlevels <- .getXlevels(Terms, m)
  fit$contrasts <- attr(X, "contrasts")
  class(fit) <- "lmmbygls"
  return(fit)
}






