#' @export
scan.EM <- function(genomecache, data, formula, model="additive",
                    chr="all", 
                    seed=1, oracle=FALSE,
                    print.locus.finished=TRUE, convergence.limit=0.0001, step.limit=1000,
                    do.augment=FALSE, just.these.loci=NULL, print.locus.fit=FALSE, ...){
 
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model="additive"))
  
  data.and.K <- make.processed.data(formula=formula, data=data, cache.subjects=cache.subjects, K=NULL, impute.on="SUBJECT.NAME")
  data <- data.and.K$data

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
  old.data <- data
  
  ###### Augmentation
  if(do.augment){
    augment.n <- ifelse(model=="additive", num.founders, num.founders + choose(num.founders, 2))
    data <- make.simple.augment.data(data=data, augment.n=augment.n)
  }
  
  ###### Null model
  ## No kinship effect - weights or no weights
  fit0 <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, use.par="h2", fix.par=0, weights=NULL, brute=FALSE)
  
  LOD.vec <- p.vec <- df <- rep(NA, length(loci))
  null.data <- data
  
  for(i in 1:length(loci)){
    X.add <- h$getLocusMatrix(loci[i], model="additive", subjects=data$SUBJECT.NAME[1:original.n])
    colnames(X.add) <- gsub(pattern="/", replacement=".", x=colnames(X.add), fixed=TRUE)
    locus.formula <- make.EM.alt.formula(formula=formula, X=X.add)
    
    X.full <- h$getLocusMatrix(loci[i], model="full", subjects=data$SUBJECT.NAME[1:original.n])
    colnames(X.full) <- gsub(pattern="/", replacement=".", x=colnames(X.full), fixed=TRUE)
    
    if(do.augment){
      X.names <- rownames(X.add)
      if(model=="additive"){
        X.add <- rbind(X.add, 2*diag(augment.n))
      }
      if(model=="full"){
        X.add <- rbind(X.add, diag(augment.n))
      }
      rownames(X.add) <- c(X.names, paste0("augment.obs", 1:augment.n))
    }
      
    ROP.data <- cbind(null.data, X.add)
    ## Using ROP to get starting values
    ROP.results <- lmmbygls(formula=locus.formula, data=ROP.data, 
                            eigen.K=NULL, K=NULL,
                            use.par="h2", fix.par=0)
    coefficients <- ROP.results$coefficients
    coefficients[is.na(coefficients)] <- 0
    sigma2 <- ROP.results$sigma2.mle
    
    expanded.diplotype.data <- get.expanded.X(X.full, coefficients=NULL)
    expanded.dosage.data <- rotate.full.to.add.data(diplotype.data=expanded.diplotype.data)
    X <- as.matrix(expanded.dosage.data[,-1])
    colnames(X) <- names(expanded.dosage.data)[-1]
    rownames(X) <- expanded.dosage.data[,1]
    X <- apply.contrast(X)
    expanded.y <- get.expanded.y(y=data[,all.vars(formula)[1]], num.groups=36)
    start.weights <- get.weights.from.diplotype.prob.data(X.full)
    
    steps <- 0
    has.converged <- FALSE
    max.dif.vec <- converge.vec <- NULL
    set.seed(seed)
    while(!has.converged){
      previous.coefficients <- coefficients
      previous.sigma2 <- sigma2
      weights <- E.step(start.weights=start.weights, y=expanded.y, X=X, coefficients=coefficients, sigma2=sigma2)
      WLS.fit <- M.step(weights=weights, y=expanded.y, X=X)
      coefficients <- WLS.fit$coefficients
      sigma2 <- WLS.fit$sigma2
      
      coefficient.dif <- get.coef.dif(coefficients=coefficients, previous.coefficients=previous.coefficients)
      max.dif <- max(abs(c(coefficient.dif, sigma2 - previous.sigma2)))
      max.dif.vec <- c(max.dif.vec, max.dif)
      if(max.dif < convergence.limit){ 
        has.converged <- TRUE 
        converge.vec <- c(converge.vec, TRUE)
      }
      steps <- steps + 1
      
      cat(paste0("step: ", steps, ", max diff: ", max.dif, "\n"))
      if(steps == step.limit){
        has.converged <- TRUE
        converge.vec <- c(converge.vec, FALSE)
      }
    }
    alt.logLik <- get.mixture.likelihood(start.weights=start.weights, expanded.X=X, 
                                         expanded.y=expanded.y, coefficients=coefficients, sigma2=sigma2)
    im.rank <- sum(!is.na(coefficients))
    
    df[i] <- im.rank
    LOD.vec[i] <- log10(exp(alt.logLik - fit0$logLik))
    p.vec[i] <- pvalue.per.locus.im(im.logLik=alt.logLik, im.rank=im.rank, fit0=fit0)
    if(print.locus.fit){ cat(paste("locus", i, "out of", length(loci)), "\n") }
  }
  names(LOD.vec) <- names(p.vec) <- names(df) <- loci
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 MI.LOD=NULL,
                 MI.p.value=NULL,
                 df=df,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 fit0.REML=NULL,
                 y=data$y,
                 formula=formula.string,
                 model.type=model)
  return(output)
}


