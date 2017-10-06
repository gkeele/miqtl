
#' @export
condition.out.locus.for.scan <- function(locus, new.outcome.name="new.y",
                                         genomecache, data, 
                                         formula, K=NULL,
                                         model=c("additive", "full"), locus.as.fixed=TRUE,
                                         use.par="h2", use.multi.impute=TRUE, num.imp=11, brute=TRUE, use.fix.par=TRUE, 
                                         seed=1, pheno.id="SUBJECT.NAME", geno.id="SUBJECT.NAME",
                                         weights=NULL, do.augment=FALSE, use.full.null=FALSE, added.data.points=1, 
                                         ...){
  model <- model[1]

  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  locus.chr <- h$getChromOfLocus(locus)

  ## Case where there are replicates and K is not specified, 
  ## then K is forced to be identity
  if(pheno.id != geno.id & is.null(K)){
    K <- diag(unique(data[,geno.id]))
    rownames(K) <- colnames(K) <- unique(data[,geno.id])
  }
  
  cache.subjects <- rownames(h$getLocusMatrix(h$getLoci()[1], model="additive"))
  data.and.K <- make.processed.data(formula=formula, data=data, 
                                    cache.subjects=cache.subjects, K=K, 
                                    pheno.id=pheno.id, geno.id=geno.id)
  data <- data.and.K$data
  K <- data.and.K$K
  ## Check that some IDs are matching between data and K
  if(!is.null(K)){
    if(all(dim(K) == 0)){
      stop("The colnames and rownames of K do not match the specified ID geno.id in the data", call.=FALSE)
    }
  }
  
  cache.subjects <- unique(as.character(data[,geno.id]))
  if(!is.null(weights)){ weights <- weights[as.character(data[,pheno.id])] }
  
  augment.indicator <- NULL
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
  original.n <- nrow(data)
  old.data <- data
  
  ###### Augmentation
  if(do.augment){
    augment.n <- ifelse(model=="additive", num.founders, num.founders + choose(num.founders, 2))
    augment.indicator <- c(rep(0, original.n), rep(1, augment.n))
    if(!use.full.null){
      data <- make.simple.augment.data(data=data, formula=formula, augment.n=augment.n)
      data <- data.frame(data, augment.indicator=augment.indicator)
      K <- make.simple.augment.K(K=K, augment.n=augment.n)
    }
    if(use.full.null){
      no.augment.K <- K
      K <- make.full.null.augment.K(K=no.augment.K, original.n=original.n, augment.n=augment.n)
      data <- make.full.null.augment.data(formula=formula, data=data, no.augment.K=no.augment.K, use.par=use.par, brute=brute,
                                          original.n=original.n, augment.n=augment.n, weights=weights)
    }
    weights <- make.augment.weights(data=data, weights=weights, augment.n=augment.n, added.data.points=added.data.points)
  }
  ###### Null model
  ## check for LMER notation
  use.lmer <- check.for.lmer.formula(null.formula)
  
  ####################### CHECKS
  ###### check that full directory has data
  if(model == "full" | use.multi.impute){
    if(!file.exists(paste0(genomecache, "/full/chr", locus.chr, "/data/", locus, ".RData"))){
      stop("Full model probabilities not available in genome cache, only additive ROP can be fit", call.=FALSE)
    }
  }
  ###### check that both sparse random effect and kinship random effect are not being specified together
  if(use.lmer & !is.null(K)){
    stop("Cannot use LMER sparse random effects AND a non-sparse random effect", call.=FALSE)
  }
  ###### check that both lmer and random locus effect are not being specified together
  if(use.lmer & !locus.as.fixed){
    stop("Cannot use LMER sparse random effects AND fit locus effect as random", call.=FALSE)
  }
  ###### check that use.fix.par=TRUE if locus.as.fixed=FALSE
  if(!use.fix.par & !locus.as.fixed){
    cat("standard ANOVA F-test not valid with mixed effect model, switching to LRT\n")
    use.fix.par <- TRUE
  }
  ###### Null model fits
  if(use.lmer){
    fit0 <- lmmbylmer(formula=null.formula, data=data, REML=FALSE, weights=weights)
    fix.par <- NULL
  }
  else{
    ## No kinship effect - weights or no weights
    if(is.null(K)){
      fit0 <- lmmbygls(formula=null.formula, data=data, eigen.K=NULL, K=NULL, pheno.id=pheno.id,
                       use.par="h2", fix.par=0, weights=weights, brute=brute)
    }
    ## Kinship effect - weights or no weights
    else{
      ###### Handling replicates
      if(pheno.id != geno.id){
        Z <- model.matrix(process.random.formula(geno.id=geno.id), data=data)
        K <- Z %*% K %*% t(Z)
        rownames(K) <- colnames(K) <- as.character(data[,pheno.id])
      }
      ###### Handling constant weights
      if(!is.null(weights)){
        J <- weights^(1/2) * t(weights^(1/2) * K)
        eigen.J <- process_eigen_decomposition(eigen.decomp=eigen(J))
        fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.J, K=K, use.par=use.par, weights=weights, brute=brute)
      }
      else{
        eigen.K <- process_eigen_decomposition(eigen.decomp=eigen(K))
        fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.K, K=K, use.par=use.par, weights=weights, brute=brute)
      }
    }
    ####### EMMA or EMMAX  
    if(use.fix.par){
      fix.par <- ifelse(locus.as.fixed, fit0$h2, fit0.REML$h2)
    }
    if(!use.fix.par){
      fix.par <- NULL
    }
  }
  null.data <- data
  
  ## Prepping link between phenotype and genotype (necessary for imputation in multiple imputations)
  impute.map <- data.frame(data[,pheno.id], data[,geno.id])
  names(impute.map) <- c(pheno.id, geno.id)
  non.augment.subjects <- as.character(data[,geno.id])[grep(pattern="augment", x=as.character(data[,geno.id]), invert=TRUE)]
  
  ## More efficient - does not require the formula/data processing steps
  y <- data$y
  
  ## Alternative model
  if(use.multi.impute){
    diplotype.prob.matrix <- h$getLocusMatrix(locus, model="full", subjects=non.augment.subjects)
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
      sample.names <- rownames(diplotype.prob.matrix)
      diplotype.prob.matrix <- rbind(diplotype.prob.matrix, augment.matrix)
      rownames(diplotype.prob.matrix) <- c(sample.names, paste0("augment.obs", 1:augment.n))
    }
    if(locus.as.fixed){ fit0.for.mi <- fit0 }
    else{
      stop("Currently cannot take residuals from random effect fitting of QTL", call.=FALSE)
      #fit0.for.mi <- fit0.REML 
    }
    fit1 <- multi.imput.lmmbygls(formula=formula,
                                 y=y, X.probs=diplotype.prob.matrix,
                                 weights=weights, locus.as.fixed=locus.as.fixed,
                                 model=model, founders=founders, pheno.id=pheno.id, num.imp=num.imp,
                                 use.lmer=use.lmer, impute.map=impute.map,
                                 use.par=use.par, fix.par=fix.par, fit0=fit0.for.mi, do.augment=do.augment, 
                                 brute=brute, seed=seed) 
  }
  else{ ## ROP
    X <- h$getLocusMatrix(locus, model=model, subjects=non.augment.subjects)
    keep.col <- 1:ncol(X)
    if(locus.as.fixed){
      max.column <- which.max(colSums(X, na.rm=TRUE))[1]
      keep.col <- keep.col[keep.col != max.column]
      X <- X[,keep.col]
    }
    colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
    locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
    if(do.augment){
      X.names <- rownames(X)
      if(model=="additive"){
        X <- rbind(X, 2*diag(augment.n)[,keep.col])
      }
      if(model=="full"){
        X <- rbind(X, diag(augment.n)[,keep.col])
      }
      rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
    }
    if(use.lmer){
      data <- cbind(null.data, X)
      fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
    }
    else{
      if(locus.as.fixed){
        X <- cbind(fit0$x, X)
        fit1 <- lmmbygls(formula=locus.formula, 
                         y=y, X=X,
                         eigen.K=fit0$eigen.K, K=fit0$K, weights=weights,
                         use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                         brute=brute)
        qtl.predictor <- regress.out.qtl(fit1=fit1, null.formula=null.formula, alt.formula=locus.formula,
                                         locus.as.fixed=locus.as.fixed)
      }
      else{
        fit1 <- lmmbygls.random(formula=null.formula, pheno.id=pheno.id,
                                y=y, X=fit0$x,
                                eigen.K=fit0$eigen.K, K=fit0$K, Z=X, weights=weights,
                                use.par="h2", null.h2=fix.par,
                                brute=brute)
      }
    }
    new.y <- y - qtl.predictor
    return.data <- null.data
    return.data[,new.outcome.name] <- new.y
  }
  return(return.data)
}


regress.out.qtl <- function(fit1, null.formula, alt.formula, null.data,
                            locus.as.fixed){
  qtl.vars <- all.vars(alt.formula)[!(all.vars(alt.formula) %in% all.vars(null.formula))]
  if(locus.as.fixed){
    qtl.predictor <- fit1$x[,qtl.vars] %*% fit1$coefficients[qtl.vars]
  }
  return(qtl.predictor)
}
