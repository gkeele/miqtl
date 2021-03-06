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

#' Run a haplotype-based genome scan from probabilities stored in a genome cache directory
#'
#' This function primarily takes a formula, data frame, and genome cache to run a genome scan.
#'
#' @param formula The lm/lmer-style formula for the model to fit. Variables must correspond to columns in the data.frame specified in the data
#' argument. If no data.frame is provided, and rather a y vector and X matrix, the formula is still passed to the resulting output fit object.
#' @param data DEFAULT: NULL. The data.frame that contains the variables included in the model. If no data.frame is specified, the expectation
#' is that a further process y vector and X matrix is provided. This allows for greater computational efficiency, side-stepping the need to pull
#' the quantities from the data.frame based on the formula.
#' @param y DEFAULT: NULL. y is an outcome vector that can be specified instead of a data.frame.
#' @param X DEFAULT: NULL. X is a design matrix that can be specified instead of a data.frame.
#' @param K DEFAULT: NULL. K is the covariance matrix, commonly a realized genetic relationship matrix or kinship matrix. NULL is interpreted as
#' no variance component should be fit. K should have row and column names that match the pheno.id column in the data.frame or the order of y and X.
#' @param eigen.K DEFAULT: NULL. Is the eigendecomposition of K. If NULL and K is non-null, it will be computed. Pre-computing in something like a
#' genome scan (as in scan.h2lmm()) can save computation time.
#' @param fix.par DEFAULT: NULL. This can be the ML or REML estimate of h2 from a null model. Using this within a genome scan results in an EMMAX-like
#' scan, which is the standard. It saves time because the parameter estimates will be analytical now rather than requiring an optimization step. 
#' @param M DEFAULT: NULL. The GLS multiplier matrix. Useful in the context of genome scans with EMMAX-like specifications, as it saves computation time.
#' @param logDetV DEFAULT: NULL. Also useful within the context of a genome scan with EMMAX.
#' @param weights DEFAULT: NULL. Allows for the unstructured error to have diagonal but non-identity covariance structure. More intuitively, it allows 
#' observations to be differentially weighted. The default of NULL is equivalent to a vector of ones. Vector should have names that match pheno.id column
#' in data.frame or match the order of y and X.
#' @param pheno.id DEFAULT: "SUBJECT.NAME". The column in data that corresponds to the individual ID.
#' @param use.par DEFAULT: "h2". Specifies that the h2 should be optimized with respect to the likelihood (ML). Alternatively, "h2.REML" can be specified,
#' which optimizes h2 in terms of the residual likelihood (REML).
#' @param brute DEFAULT: TRUE. Internally optim() is used to find ML or REML estimates. It does not explicitly check the boundaries of h2. This option forces
#' it to check h2=0, which can occur.
#' @param subset Option to subset the data.
#' @param na.action Determines how NAs are handled.
#' @param method DEFAULT: "qr". Option used by lm.fit after GLS process is completed.
#' @param model DEFAULT: TRUE. Return model output. Can be turned off for efficiency within a scan.
#' @param contrasts DEFAULT: NULL.
#' @param verbose DEFAULT: FALSE.
#' @export
#' @examples lmmbygls()
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
        eigen.K <- eigen(K, symmetric=TRUE)
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
    
    logDetXtX <- log(det(crossprod(X)))
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
  big.K <- Z %*% tcrossprod(use.K, Z)
  rownames(big.K) <- colnames(big.K) <- as.character(use.data[,pheno.id])
  lmmbygls.fit <- lmmbygls(formula=formula, data=use.data, K=big.K,
                           pheno.id=pheno.id, ...)
  return(lmmbygls.fit)
}



