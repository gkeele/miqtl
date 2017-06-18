lmmbygls.random <- function(formula, data, Z, 
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
  if(method == "model.frame"){
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
  h2.fit <- function(h2, logLik.only=TRUE, ...){
    H <- K*h2 + diag(d*null.h2*(1 - h2) + 1 - null.h2*(1 - h2) - h2)
    chol.H <- chol(H)
    M <- t(solve(chol.H))
    fit <- gls.fit(X=X, y=y, M=M, ...)
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






