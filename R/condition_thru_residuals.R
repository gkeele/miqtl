#' Add the residuals to a data.frame after regressing out a specified QTL effect
#'
#' This function sets up conditional scans via regressing out QTL effects. The output is a new data.frame
#' with a residual outcome to be used as phenotype in subsequent genome scans.
#'
#' @param locus The locus for which the QTL effect will be regressed out. The locus must be contained in the
#' genomecache.
#' @param new.outcome.name DEFAULT: "new.y". The name of the new residual outcome variable to be included in the output data set.
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often with the individual-level ID named "SUBJECT.NAME", though others
#' can be specified with pheno.id.
#' @param formula An lm style formula with functions of outcome and covariates contained in data frame.
#' @param K DEFAULT: NULL. A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
#' based on SNP genotypes or the founder haplotype probabilities. Colnames and rownames should match
#' the SUBJECT.NAME column in the data frame. If no K matrix is specified, either lmer is used (if sparse random effects
#' are included in the formula) or a fixed effect model (equivalent to lm).
#' @param model DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities.
#' @param locus.as.fixed DEFAULT: TRUE. If TRUE, the locus effect is fit as fixed effect. If FALSE, it is
#' fit as a random effect.
#' @param use.par DEFAULT: "h2". The parameterization of the likelihood to be used. 
#' @param use.multi.impute DEFAULT: TRUE. This option specifies whether to use ROP or multiple imputations.
#' @param num.imp DEFAULT: 11. IF multiple imputations are used, this specifies the number of imputations to perform.
#' @param brute DEFAULT: TRUE. During the optimization to find maximum likelihood parameters, this specifies checking the
#' boundary of h2=0. Slightly less efficient, but otherwise the optimization procedure will not directly check
#' the boundary.
#' @param use.fix.par DEFAULT: TRUE. This specifies an approximate fitting of mixed effect model (Kang et al. 2009). Much
#' more efficient, as the optimization of h2 only needs to be performed once for the null model rather than every locus. 
#' Technically less powerful, though in practice it has proven to be almost equal to the exact procedure.
#' @param seed DEFAULT: 1. Multiple imputations involve a sampling process of the diplotypes, thus a seed is necessary
#' to produce the same results over multiple runs and different machines.
#' @param pheno.id DEFAULT: "SUBJECT.NAME". The is the individual-level ID that is associated with data points in the phenotype
#' data. Generally this should be unique for each data point.
#' @param geno.id DEFAULT: "SUBJECT.NAME". The default represents the situation that each genome is unique. Specifying some other
#' column allows for replicate genomes, such as in the CC or CC-RIX.
#' @param weights DEFAULT: NULL. If unspecified, individuals are equally weighted. This option allows for a weighted analysis 
#' when using the mean of multiple individuals with the same genome.
#' @param do.augment DEFAULT: FALSE. Augments the data with null observations for genotype groups. This is an approximately Bayesian 
#' approach to applying a prior to the data, and can help control highly influential data points.
#' @param use.full.null DEFAULT: FALSE. Draws augmented data points from the null model. This allows for the inclusion of null data points
#' that do not influence the estimation of other model parameters as much.
#' @param added.data.points DEFAULT: 1. If augment weights are being used, this specifies how many data points should be added in total.
#' @export
#' @examples condition.out.locus.for.scan()
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
    fit1 <- multi.imput.lmmbygls(formula=formula, null.formula=null.formula,
                                 y=y, X.probs=diplotype.prob.matrix,
                                 weights=weights, locus.as.fixed=locus.as.fixed,
                                 model=model, founders=founders, pheno.id=pheno.id, num.imp=num.imp,
                                 use.lmer=use.lmer, impute.map=impute.map,
                                 use.par=use.par, fix.par=fix.par, fit0=fit0.for.mi, do.augment=do.augment, 
                                 brute=brute, seed=seed, return.qtl.predictor=TRUE)
    qtl.predictor <- rowMeans(fit1$qtl.predictor)
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
  }
  new.y <- y - qtl.predictor
  return.data <- null.data
  return.data[,new.outcome.name] <- new.y
  return(return.data)
}


regress.out.qtl <- function(fit1, null.formula, alt.formula,
                            locus.as.fixed){
  qtl.vars <- all.vars(alt.formula)[!(all.vars(alt.formula) %in% all.vars(null.formula))]
  if(locus.as.fixed){
    qtl.vars <- names(fit1$coefficients[qtl.vars][!is.na(fit1$coefficients[qtl.vars])])
    qtl.predictor <- fit1$x[,qtl.vars] %*% fit1$coefficients[qtl.vars]
  }
  return(qtl.predictor)
}
