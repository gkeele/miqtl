#' Run a haplotype-based genome scan from probabilities stored in a genome cache directory
#'
#' This function primarily takes a formula, data frame, and genome cache to run a genome scan.
#'
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
#' @param p.value.method DEFAULT: "LRT". "LRT" specifies a likelihood ratio test, which is flexible to testing fixed
#' effects in fixed and mixed effect models. "ANOVA" specifies an F-test, which is only valid in fixed effect models.
#' ANOVA is more conservative in models with low sample sizes, where the asymptotic theory underlying the LRT does not
#' hold.
#' @param use.par DEFAULT: "h2". The parameterization of the likelihood to be used. 
#' @param use.multi.impute DEFAULT: TRUE. This option specifies whether to use ROP or multiple imputations.
#' @param num.imp DEFAULT: 11. IF multiple imputations are used, this specifies the number of imputations to perform.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
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
#' @param just.these.loci DEFAULT: NULL. Specifies a reduced set of loci to fit. If loci is just one locus, the alternative model fit
#' will also be output as fit1.
#' @param print.locus.fit DEFAULT: FALSE. If TRUE, prints out how many loci have been fit currently.
#' @param debug.single.fit DEFAULT: FALSE. If TRUE, a browser() call is activated after the first locus is fit. This option
#' allows developers to more easily debug while still using the actual R package.
#' @export
#' @examples scan.h2lmm()
scan.h2lmm <- function(genomecache, data, 
                       formula, K=NULL,
                       model=c("additive", "full"),
                       p.value.method=c("LRT", "ANOVA"),
                       use.par="h2", use.multi.impute=TRUE, num.imp=11, chr="all", brute=TRUE, use.fix.par=TRUE, 
                       seed=1, pheno.id="SUBJECT.NAME", geno.id="SUBJECT.NAME",
                       weights=NULL, do.augment=FALSE, use.full.null=FALSE, added.data.points=1, 
                       just.these.loci=NULL,
                       print.locus.fit=FALSE, debug.single.fit=FALSE,
                       ...){
  model <- model[1]
  p.value.method <- p.value.method[1]
  
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model="additive"))
  data.and.K <- make.processed.data(formula=formula, data=data, 
                                    cache.subjects=cache.subjects, K=K, 
                                    pheno.id=pheno.id, geno.id=geno.id)
  data <- data.and.K$data
  K <- data.and.K$K
  cache.subjects <- unique(as.character(data[,geno.id]))
  if(!is.null(weights)){ weights <- weights[as.character(data[,pheno.id])] }

  loci.chr <- h$getChromOfLocus(loci)
  if(chr[1] != "all"){
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
  
  ###### check p-value method
  if(use.lmer & p.value.method == "ANOVA"){
    cat("ANOVA not currently supported with our implementation of LMER, swithcing to LRT\n")
    p.value.method <- "LRT"
  }
  if(p.value.method == "ANOVA" & !is.null(K)){
    cat("standard ANOVA F-test not valid with mixed effect model, swithcing to LRT\n")
    p.value.method <- "LRT"
  }
  
  ###### check that both sparse random effect and kinship random effect are not being specified together
  if(use.lmer & !is.null(K)){
    stop("Cannot use LMER sparse random effects AND a non-sparse random effect", call.=FALSE)
  }
  
  ###### Null model fits
  if(use.lmer){
    fit0 <- lmmbylmer(formula=null.formula, data=data, REML=FALSE, weights=weights)
    fit0.REML <- lmmbylmer(formula=null.formula, data=data, REML=TRUE, weights=weights)
    fix.par <- NULL
  }
  else{
    ## No kinship effect - weights or no weights
    if(is.null(K)){
      fit0 <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, pheno.id=pheno.id,
                       use.par="h2", fix.par=0, weights=weights, brute=brute)
      fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, pheno.id=pheno.id,
                            use.par="h2.REML", fix.par=0, weights=weights, brute=brute)
    }
    ## Kinship effect - weights or no weights
    else{
      ###### Handling replicates
      if(pheno.id != geno.id){
        Z <- model.matrix(process.random.formula(geno.id=geno.id), data=data)
        K <- Z %*% K %*% t(Z)
        rownames(K) <- colnames(K) <- as.character(data[,pheno.id])
      }
      ###### Handling constant weights at all loci
      if(!is.null(weights)){
        J <- weights^(1/2) * t(weights^(1/2) * K)
        eigen.J <- process_eigen_decomposition(eigen.decomp=eigen(J))
        fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.J, K=K, use.par=use.par, weights=weights, brute=brute)
        fit0.REML <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.J, K=K, use.par="h2.REML", weights=weights, brute=brute)
      }
      else{
        eigen.K <- process_eigen_decomposition(eigen.decomp=eigen(K))
        fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.K, K=K, use.par=use.par, weights=weights, brute=brute)
        fit0.REML <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=weights, brute=brute)
      }
    }
    ####### EMMA or EMMAX  
    if(use.fix.par){
      fix.par <- fit0$h2
    }
    if(!use.fix.par){
      fix.par <- NULL
    }
  }
  MI.LOD <- MI.p.value <- NULL
  LOD.vec <- p.vec <- df <- rep(NA, length(loci))
  null.data <- data
  
  ## Prepping link between phenotype and genotype (necessary for imputation in multiple imputations)
  impute.map <- data.frame(data[,pheno.id], data[,geno.id])
  names(impute.map) <- c(pheno.id, geno.id)
  non.augment.subjects <- as.character(data[,geno.id])[grep(pattern="augment", x=as.character(data[,geno.id]), invert=TRUE)]

  for(i in 1:length(loci)){
    if(use.multi.impute){
      if(i == 1){ # only at the beginning
        MI.LOD <- MI.p.value <- matrix(NA, nrow=num.imp, ncol=length(loci))
      }
      diplotype.prob.matrix <- h$getLocusMatrix(loci[i], model="full", subjects=non.augment.subjects)
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
      fit1 <- multi.imput.lmmbygls(num.imp=num.imp, data=data, formula=formula, founders=founders,
                                   diplotype.probs=diplotype.prob.matrix, pheno.id=pheno.id,
                                   model=model, p.value.method=p.value.method, 
                                   use.lmer=use.lmer, impute.map=impute.map,
                                   use.par=use.par, fix.par=fix.par, fit0=fit0, do.augment=do.augment, 
                                   brute=brute, seed=seed, weights=weights) 
      MI.LOD[,i] <- fit1$LOD
      MI.p.value[,i] <- fit1$p.value
      LOD.vec[i] <- median(fit1$LOD)
      p.vec[i] <- median(fit1$p.value)
    }
    if(!use.multi.impute){
      X <- h$getLocusMatrix(loci[i], model=model, subjects=non.augment.subjects)
      max.column <- which.max(colSums(X, na.rm=TRUE))[1]
      X <- X[,-max.column]
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
      locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
      if(do.augment){
        X.names <- rownames(X)
        if(model=="additive"){
          X <- rbind(X, 2*diag(augment.n)[,-max.column])
        }
        if(model=="full"){
          X <- rbind(X, diag(augment.n)[,-max.column])
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
      }
      
      data <- cbind(null.data, X)
      if(use.lmer){
        fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
        LOD.vec[i] <- log10(exp(as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))))
        p.vec[i] <- pchisq(q=-2*(as.numeric(logLik(fit0)) - as.numeric(logLik(fit1))), df=length(fixef(fit1))-length(fixef(fit0)), lower.tail=FALSE)
      }
      else{
        fit1 <- lmmbygls(formula=locus.formula, data=data, pheno.id=pheno.id,
                         eigen.K=fit0$eigen.K, K=fit0$K, 
                         use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                         brute=brute, 
                         weights=weights)
        LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
        p.vec[i] <- get.p.value(fit0=fit0, fit1=fit1, method=p.value.method)
        df[i] <- fit1$rank
      }
    }
    if(debug.single.fit){ browser() }
    if(print.locus.fit){ cat(paste("locus", i, "out of", length(loci)), "\n") }
  }
  names(LOD.vec) <- names(p.vec) <- names(df) <- loci
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 MI.LOD=MI.LOD,
                 MI.p.value=MI.p.value,
                 df=df,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 fit0.REML=fit0.REML,
                 y=data$y,
                 formula=formula.string,
                 model.type=model,
                 p.value.method=p.value.method,
                 impute.map=impute.map)
  if(length(just.these.loci) == 1){ output$fit1 <- fit1 }
  if(pheno.id != geno.id & !is.null(K)){ rownames(Z) <- as.character(data[, pheno.id]); output$Z <- Z }
  return(output)
}