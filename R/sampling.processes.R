#' Returns a matrix of outcome samples, either permutations or from the null model of no locus effect
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of outcome samples, either permutations or
#' from the null model of no locus effect.
#'
#' @param scan.object A scan.h2lmm() object.
#' @param model.type DEFAULT: "null". "null" specifies sampling processes from the null model. "alt" specifies sampling processes
#' from the alternative model.
#' @param method DEFAULT: "bootstrap". "bootstrap" specifies parametric bootstraps from the given model. "permutation" specifies
#' parametric permutations that can respect the structure of the data. Permutations are more appropriate if the data have highly
#' influential data points.
#' @param use.REML DEFAULT: TRUE. Determines whether the variance components for the parametric sampling are 
#' based on maximizing the likelihood (ML) or the residual likelihood (REML).
#' @param use.BLUP DEFAULT: FALSE.This results in the BLUP value of the polgyene effect (assuming a GRM has been given) is used,
#' rather than sampled. This reduces the variation seen across sampling, which can result in narrower positional confidence 
#' intervals.
#' @param num.samples The number of parametric bootstrap samples to return.
#' @param seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples generate.sample.outcomes.matrix()
generate.sample.outcomes.matrix <- function(scan.object, model.type=c("null", "alt"), 
                                            method=c("bootstrap", "permutation"), use.REML=TRUE, 
                                            use.BLUP=FALSE, num.samples, seed=1){
  model.type <- model.type[1]
  method <- method[1]
  
  if(model.type == "null"){ fit <- scan.object$fit0; locus <- NULL }
  if(model.type == "alt"){ fit <- scan.object$fit1; locus <- scan.object$loci }
  fit0.REML <- scan.object$fit0.REML
  if(class(fit) != "lmerMod"){
    Xb <- fit$x %*% fit$coefficients
    n <- nrow(fit$x)
    K <- fit$K
    weights <- fit$weights
    return.weights <- weights
    if(is.null(weights)){ 
      weights <- rep(1, nrow(K)) 
    }
    if(use.REML){
      if(is.null(K)){
        tau2 <- 0
        sigma2 <- fit$sigma2.mle*(n/(n - 1))
      }
      else{
        tau2 <- fit0.REML$tau2.mle
        sigma2 <- fit0.REML$sigma2.mle
      }
    }
    else{
      tau2 <- fit$tau2.mle
      sigma2 <- fit$sigma2.mle  
    }
    sim.y.matrix <- matrix(NA, nrow=n, ncol=num.samples)
    
    if(is.null(K)){
      u <- rep(0, n)
    }
    else{
      original.K <- K
      impute.map <- scan.object$impute.map
      K <- reduce.large.K(large.K=K, impute.map=impute.map)
      if(use.BLUP){
        X <- fit$x
        Sigma <- original.K*tau2 + diag(1/weights)*sigma2
        inv.Sigma <- solve(Sigma)
        u.BLUP <- (original.K*tau2) %*% inv.Sigma %*% (diag(nrow(original.K)) - X %*% solve(t(X) %*% inv.Sigma %*% X) %*% t(X) %*% inv.Sigma) %*% fit$y  
      }
    }
    
    set.seed(seed)
    for(i in 1:num.samples){
      if(!is.null(K)){
        ## Handling potential replicates
        if(use.BLUP){
          u <- u.BLUP
        }
        else{
          u <- c(mnormt::rmnorm(1, mean=rep(0, nrow(K)), varcov=K*tau2))
        }
        names(u) <- unique(impute.map[,2])
        u <- u[impute.map[,2]]
      }
      if(is.null(weights)){
        e <- rnorm(n=n, mean=0, sd=sqrt(sigma2))
      }
      else{
        e <- c(mnormt::rmnorm(1, mean=rep(0, n), varcov=diag(1/weights)*sigma2))
      }
      y.sample <- Xb + u + e
      if(method == "bootstrap"){
        sim.y.matrix[,i] <- y.sample
      }
      if(method == "permutation"){
        perm.y.ranks <- order(y.sample)
        sim.y.matrix[,i] <- fit$y[perm.y.ranks]
      }
    }
    rownames(sim.y.matrix) <- names(fit$y)
  }
  else{
    stop("Need to add lmer-based functionality!!")
  }
  sim.threshold.object <- list(y.matrix=sim.y.matrix,
                               formula=scan.object$formula,
                               weights=return.weights,
                               K=K,
                               method=method,
                               impute.map=scan.object$impute.map,
                               locus=locus)
  return(sim.threshold.object)
}

### Support function that can take a large K with replicates rows/columns, and reduce it down
reduce.large.K <- function(large.K, impute.map){
  map.order <- match(impute.map[,1], table=colnames(large.K))
  impute.map <- impute.map[map.order,]
  colnames(large.K) <- rownames(large.K) <- impute.map[,2]
  K <- large.K[unique(as.character(impute.map[,2])), unique(as.character(impute.map[,2]))]
  return(K)
}

#' Runs threshold scans from a matrix of outcomes, either parametric bootstraps from the null model or permutations
#'
#' This function takes an object produced from either generate.null.bootstrap.matrix() or generate.perm.matrix(), and 
#' runs genome scans on the outcomes contained in them.
#'
#' @param sim.threshold.object An object created by either generate.null.bootstrap.matrix() or generate.perm.matrix().
#' @param keep.full.scans DEFAULT: TRUE. Returns full genome scans for every outcome sample in the sim.threshold.object. Can be used
#' for visualization of the procedure, but greatly increases the size of the output object.
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often the individual-level ID named "SUBJECT.NAME".
#' @param model DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities.
#' @param use.multi.impute DEFAULT: TRUE. This option specifies whether to use ROP or multiple imputations.
#' @param num.imp DEFAULT: 11. IF multiple imputations are used, this specifies the number of imputations to perform.
#' @param chr DEFAULT: "all". The chromosomes to conduct scans over.
#' @param just.these.loci DEFAULT: NULL. Specifies a reduced set of loci to fit. If loci is just one locus, the alternative model fit
#' will also be output as fit1.
#' @param scan.seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples run.threshold.scans()
run.threshold.scans <- function(sim.threshold.object, keep.full.scans=TRUE,
                                genomecache, data,
                                model=c("additive", "full"),
                                use.multi.impute=TRUE, num.imp=11, chr="all", just.these.loci=NULL, 
                                scan.seed=1, ...){
  y.matrix <- sim.threshold.object$y.matrix
  formula <- sim.threshold.object$formula
  weights <- sim.threshold.object$weights
  K <- sim.threshold.object$K
  pheno.id <- names(sim.threshold.object$impute.map)[1]
  geno.id <- names(sim.threshold.object$impute.map)[2]
  
  num.scans <- ncol(y.matrix)
  
  h <- DiploprobReader$new(genomecache)
  loci <- h$getLoci()
  loci.chr <- h$getChromOfLocus(loci)
  if(chr[1] != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  full.p <- full.lod <- these.chr <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- full.lod <- matrix(NA, nrow=num.scans, ncol=length(loci))
    colnames(full.p) <- colnames(full.lod) <- loci
    these.chr <- h$getChromOfLocus(loci)
    these.pos <- list(Mb=h$getLocusStart(loci, scale="Mb"),
                      cM=h$getLocusStart(loci, scale="cM"))
  }
  min.p <- max.lod <- rep(NA, num.scans)
  
  iteration.formula <- formula(paste0("new_y ~ ", unlist(strsplit(formula, split="~"))[-1]))
  for(i in 1:num.scans){
    new.y <- data.frame(y.matrix[,i], rownames(y.matrix))
    names(new.y) <- c("new_y", pheno.id)
    this.data <- merge(x=new.y, y=data, by=pheno.id, all.x=TRUE)
    
    this.scan <- scan.h2lmm(genomecache=genomecache, data=this.data, 
                            formula=iteration.formula, K=K, model=model,
                            use.multi.impute=use.multi.impute, num.imp=num.imp, 
                            pheno.id=pheno.id, geno.id=geno.id, seed=scan.seed,
                            weights=weights, chr=chr,
                            ...)
    if(keep.full.scans){
      full.p[i,] <- this.scan$p.value
      full.lod[i,] <- this.scan$LOD
    }
    min.p[i] <-  min(this.scan$p.value)
    max.lod[i] <- max(this.scan$LOD)
    cat("threshold scan:", i, "\n")
  }
  return(list(full.results=list(LOD=full.lod,
                                p.value=full.p,
                                chr=these.chr, 
                                pos=these.pos), 
              max.statistics=list(LOD=max.lod,
                                  p.value=min.p)))
}

################# QTL CI methods


#' Run single chromosome scans on parametric bootstrap samples from the alternative model of a particular locus
#'
#' This function runs single chromosome scans of parametric bootstrap samples from the alternative model of a particular locus. These
#' association scans can then be used to estimate a confidence interval on QTL position.
#'
#' @param sim.object A sample object. This is primarily a matrix of parametric bootstrap
#' samples from a locus.
#' @param keep.full.scans DEFAULT: TRUE. If TRUE, the full scans from each sample are kept. This allows for them
#' to be plotted out, but also increases the memory or storage needed.
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often the individual-level ID named "SUBJECT.NAME".
#' @param model DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities.
#' @param use.par DEFAULT: "h2". The parameterization of the likelihood to be used.
#' @param use.multi.impute DEFAULT: TRUE. If TRUE, use multiple imputations of genetic data. If FALSE, use ROP.
#' @param num.imp DEFAULT: 11. If multiple imputations are used, this specifies the number of imputations to perform.
#' @param brute DEFAULT: TRUE. During the optimization to find maximum likelihood parameters, this specifies checking the
#' boundaries of h2=0 and h2=1. Slightly less efficient, but otherwise the optimization procedure will not directly check
#' these values.
#' @param use.fix.par DEFAULT: TRUE. This specifies an approximate fitting of mixed effect model (Kang et al. 2009). Much
#' more efficient, as the optimization of h2 only needs to be performed once for the null model rather than every locus. 
#' Technically less powerful, though in practice it has proven to be almost equal to the exact procedure.
#' @param scan.seed DEFAULT: 1. If imputations are used, the sampling procedure is a random process. Specifying a seed allows 
#' consistent results across runs and machines. This is reset over scans, so that the same imputations are used for each scan.
#' @param scan.seed DEFAULT: 1.
#' @param do.augment DEFAULT: FALSE. Augments the data with null observations for genotype groups. This is an approximately Bayesian 
#' approach to applying a prior to the data, and can help control highly influential data points.
#' @param use.augment.weights DEFAULT: FALSE. Specify non-equal weights on the augmented data points. This allows for the inclusion of
#' augmented data points to all genotype classes while reducing their overall contribution to the data.
#' @param use.full.null DEFAULT: FALSE. Draws augmented data points from the null model. This allows for the inclusion of null data points
#' that do not influence the estimation of other model parameters as much.
#' @param added.data.points DEFAULT: 1. If augment weights are being used, this specifies how many data points should be added in total.
#' @export
#' @examples run.positional.scans()
run.positional.scans <- function(sim.object, keep.full.scans=TRUE,
                                 genomecache, data,
                                 model=c("additive", "full"),
                                 use.par="h2", use.multi.impute=TRUE, num.imp=11, brute=TRUE, use.fix.par=FALSE, 
                                 scan.seed=1, do.augment=FALSE,
                                 use.augment.weights=FALSE, use.full.null=FALSE, added.data.points=1,
                                 ...){
  model <- model[1]
  
  y.matrix <- sim.object$y.matrix
  formula <- sim.object$formula
  weights <- sim.object$weights
  K <- sim.object$K
  
  impute.map <- sim.object$impute.map
  pheno.id <- colnames(impute.map)[1]
  geno.id <- colnames(impute.map)[2]
  
  num.scans <- ncol(y.matrix)
  
  h <- DiploprobReader$new(genomecache)
  chr <- h$getChromOfLocus(sim.object$locus)
  loci <- h$getLoci()
  chr.of.loci <- h$getChromOfLocus(loci)
  loci <- loci[chr.of.loci == chr]
  
  full.results <- these.chr <- these.pos <- NULL
  if(keep.full.scans){
    full.results <- matrix(NA, nrow=num.scans, ncol=length(loci))
    colnames(full.results) <- loci
    these.chr <- rep(chr, length(loci))
    these.pos <- list(Mb=h$getLocusStart(loci, scale="Mb"),
                      cM=h$getLocusStart(loci, scale="cM"))
  }
  max.results <- rep(NA, num.scans)
  
  peak.loci.vec <- rep(NA, num.scans)
  iteration.formula <- formula(paste0("new.y ~ ", unlist(strsplit(formula, split="~"))[-1]))
  for(i in 1:num.scans){
    new.y <- data.frame(y.matrix[,i], row.names(y.matrix))
    names(new.y) <- c("new.y", pheno.id)
    this.data <- merge(x=new.y, y=data, by=pheno.id, all.x=TRUE)
    
    this.scan <- scan.h2lmm(genomecache=genomecache, data=this.data, formula=iteration.formula, K=K, model=model,
                            use.par=use.par, use.multi.impute=use.multi.impute, num.imp=num.imp, chr=chr, brute=brute, use.fix.par=use.fix.par, seed=scan.seed, do.augment=do.augment, 
                            weights=weights, use.augment.weights=use.augment.weights, use.full.null=use.full.null, added.data.points=added.data.points,
                            pheno.id=pheno.id, geno.id=geno.id)
    peak.index <- which.max(-log10(this.scan$p.value))
    peak.loci.vec[i] <- this.scan$loci[peak.index]
    if(keep.full.scans){
      full.results[i,] <- this.scan$p.value
    }
    cat("positional scan:", i, "\n")
  }
  peak.pos <- list(Mb=h$getLocusStart(loci=peak.loci.vec, scale="Mb"),
                   cM=h$getLocusStart(loci=peak.loci.vec, scale="cM"))
  return(list(full.results=list(p.values=full.results, chr=these.chr, pos=these.pos), 
              peak.loci=peak.loci.vec, 
              peak.loci.pos=peak.pos, 
              ci=list(Mb=quantile(peak.pos$Mb, probs=c(0.05, 0.95)),
                      cM=quantile(peak.pos$cM, probs=c(0.05, 0.95))), 
              chr=chr))
}

averaged.mi.parametric.bootstrap <- function(formula, data, genomecache, K,
                                             peak.locus, 
                                             num.av.over, num.bs.samples, num.imp,
                                             model, use.par="h2", fix.par=NULL, use.scan.fix.par=TRUE,
                                             do.augment=FALSE, seed=1, scale="cM"){
  full.to.dosages <- straineff.mapping.matrix()
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  eigen.K <- eigen(K)
  chol.K <- chol(K)
  
  P <- h$getLocusMatrix(locus=peak.locus, model="full", subjects=data$SUBJECT.NAME)
  
  chr.of.locus <- h$getChromOfLocus(loci=peak.locus)
  
  q.matrix <- matrix(NA, nrow=nrow(data), ncol=num.av.over)
  sigma2.vec <- tau2.vec <- rep(NA, num.av.over)
  
  set.seed(seed)
  cat("Fiting", num.av.over, "imputations to average:\n")
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
    locus.formula <- make.alt.formula()
    use.data <- cbind(data, X)
    fit1 <- lmmbygls(locus.formula, data=use.data, eigen.K=eigen.K, K=K, 
                     use.par=use.par, fix.par=fix.par, use.diplolasso=FALSE,
                     brute=FALSE, weights=NULL, use.chol=FALSE)
    cat(i, " ")
    keep <- !is.na(fit1$coefficients)
    q.matrix[,i] <- fit1$x[,keep] %*% fit1$coefficients[keep]
    sigma2.vec[i] <- fit1$sigma2.mle
    tau2.vec[i] <- fit1$tau2.mle
  }
  q.bar <- rowMeans(q.matrix)
  sigma2.bar <- mean(sigma2.vec)
  tau2.bar <- mean(tau2.vec)
  
  y.bs.matrix <- matrix(NA, nrow=nrow(data), ncol=num.bs.samples)
  for(i in 1:num.bs.samples){
    u <- t(chol.K) %*% rnorm(n=nrow(data), mean=0, sd=sqrt(tau2.bar))
    e <- rnorm(n=nrow(data), mean=0, sd=sqrt(sigma2.bar))
    y.bs.matrix[,i] <- q.bar + u + e
  }
  colnames(y.bs.matrix) <- paste0("y.bs.", 1:num.bs.samples)
  cat("Generated", num.bs.samples, "averaged bootstrap samples\n")
  cat("Beginning MI scans\n")
  
  new.data <- cbind(y.bs.matrix, data)
  peak.loci.vec <- rep(NA, num.bs.samples)
  for(i in 1:num.bs.samples){
    this.formula.string <- Reduce(paste, deparse(formula))
    this.formula <- formula(paste0("y.bs.", i, " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
    this.bs.scan <- scan.h2lmm(genomecache=genomecache,
                               formula=this.formula,
                               data=new.data,
                               K=K,
                               model=model,
                               use.par="h2", use.multi.impute=TRUE, num.imp=num.imp,
                               use.fix.par=use.scan.fix.par, chr=chr.of.locus)
    peak.index <- which.max(-log10(this.bs.scan$p.value))
    peak.loci.vec[i] <- this.bs.scan$loci[peak.index]
    cat("Finished scan of", i, "bootstrap sample out of", num.bs.samples, "\n")
  }
  peak.pos.vec <- h$getLocusStart(loci=peak.loci.vec, scale=scale)
  return(list(loci=peak.loci.vec, pos=peak.pos.vec, scale=scale, ci=quantile(peak.pos.vec, probs=c(0.05, 0.95)), chr=chr.of.locus))
}

double.averaged.mi.parametric.bootstrap <- function(formula, data, genomecache, K,
                                                    peak.locus,
                                                    num.av.over, num.first.bs.samples, num.second.bs.samples=1, num.qtl=1, num.imp,
                                                    model, use.par="h2", fix.par=NULL, use.scan.fix.par=TRUE,
                                                    do.augment=FALSE, seed=1, scale="cM"){
  
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/ci_support.R", local=TRUE)
  
  full.to.dosages <- straineff.mapping.matrix()
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  eigen.K <- eigen(K)
  chol.K <- chol(K)
  
  P <- h$getLocusMatrix(locus=peak.locus, model="full", subjects=data$SUBJECT.NAME)
  chr <- h$getChromOfLocus(loci=peak.locus)
  
  q.matrix <- matrix(NA, nrow=nrow(data), ncol=num.av.over)
  sigma2.vec <- tau2.vec <- rep(NA, num.av.over)
  
  set.seed(seed)
  cat("First BS: fitting", num.av.over, "imputations to average:\n")
  bs.phenotype.prefix <- "y.bs." # For the closure function
  num.bs.samples <- num.first.bs.samples # For the closure function
  y.bs.matrix <- simulate.from.average.over.imputations()
  cat("First BS: generated", num.first.bs.samples, "averaged bootstrap samples\n")
  cat("First BS: beginning scans\n")
  
  new.data <- cbind(y.bs.matrix, data)
  
  loci <- h$getLoci()
  these.loci <- loci[h$getChromOfLocus(loci) == chr]
  
  ## Setting up data structures to record scan statistics
  first.peak.loci.matrix <- matrix(NA, nrow=num.first.bs.samples, ncol=1)
  first.all.scans <- array(NA, dim=c(num.first.bs.samples, 1, length(these.loci)))
  second.peak.loci.matrix <- matrix(NA, nrow=num.first.bs.samples*num.second.bs.samples, ncol=1)
  second.all.scans <- array(NA, dim=c(num.first.bs.samples*num.second.bs.samples, 1, length(these.loci)))
  
  for(i in 1:num.bs.samples){
    this.formula.string <- Reduce(paste, deparse(formula))
    this.formula <- formula(paste0("y.bs.", i, " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
    this.first.bs.scan <- scan.h2lmm(genomecache=genomecache,
                                     formula=this.formula,
                                     data=new.data,
                                     K=K,
                                     model=model,
                                     use.par="h2", use.multi.impute=TRUE, num.imp=num.imp,
                                     use.fix.par=use.scan.fix.par, chr=chr)
    first.all.scans[i,1,] <- this.first.bs.scan$p.value
    first.bs.peak <- this.first.bs.scan$loci[which.max(-log10(this.first.bs.scan$p.value))]
    first.peak.loci.matrix[i,1] <- first.bs.peak
    
    cat("First BS: finished scan", i, "out of", num.first.bs.samples, "\n")
    
    # Second Bootstrap
    P <- h$getLocusMatrix(locus=first.bs.peak, model="full", subjects=data$SUBJECT.NAME)
    #chr.of.locus <- h$getChromOfLocus(loci=first.bs.peak)
    
    q.matrix <- matrix(NA, nrow=nrow(data), ncol=num.av.over)
    sigma2.vec <- tau2.vec <- rep(NA, num.av.over)
    
    set.seed(seed)
    cat("\tSecond BS: fitting", num.av.over, "imputations to average:\n")
    bs.phenotype.prefix <- "y.second.bs." # For the closure function
    num.bs.samples <- num.second.bs.samples # For the closure function
    second.y.bs.matrix <- simulate.from.average.over.imputations()
    cat("\tSecond BS: generated", num.second.bs.samples, "averaged bootstrap samples\n")
    cat("\tSecond BS: beginning scans\n")
    
    newer.data <- cbind(second.y.bs.matrix, data)
    for(j in 1:num.second.bs.samples){
      this.second.formula <- formula(paste0("y.second.bs.", j, " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
      this.second.bs.scan <- scan.h2lmm(genomecache=genomecache,
                                        formula=this.second.formula,
                                        data=newer.data,
                                        K=K,
                                        model=model,
                                        use.par="h2", use.multi.impute=TRUE, num.imp=num.imp,
                                        use.fix.par=use.scan.fix.par, chr=chr)
      second.all.scans[i,1,] <- this.second.bs.scan$p.value # works because I've locked num.second.bs.samples = 1
      second.bs.peak <- this.second.bs.scan$loci[which.max(-log10(this.second.bs.scan$p.value))]
      second.peak.loci.matrix[i,1] <- second.bs.peak
      cat("\tSecond BS: finished scan", j, "out of", num.second.bs.samples, "\n")
    }
  }
  first.peak.cm.matrix <- apply(first.peak.loci.matrix, 1, function(x) h$getLocusStart(loci=x, scale="cM"))
  first.peak.mb.matrix <- apply(first.peak.loci.matrix, 1, function(x) h$getLocusStart(loci=x, scale="Mb"))
  first.peak.pos.list <- list(cM=first.peak.cm.matrix, Mb=first.peak.mb.matrix)
  
  second.peak.cm.matrix <- apply(second.peak.loci.matrix, 1, function(x) h$getLocusStart(loci=x, scale="cM"))
  second.peak.mb.matrix <- apply(second.peak.loci.matrix, 1, function(x) h$getLocusStart(loci=x, scale="Mb"))
  second.peak.pos.list <- list(cM=second.peak.cm.matrix, Mb=second.peak.mb.matrix)
  
  pos <- list(cM=h$getLocusStart(these.loci, scale="cM"), Mb=h$getLocusStart(these.loci, scale="Mb"))
  
  return(list(single.bs=list(loci=these.loci,
                             pos=pos,
                             full.results=first.all.scans,
                             peak.loci=first.peak.loci.matrix,
                             peak.pos=first.peak.pos.list,
                             ci=quantile(first.peak.pos.list[[scale]], probs=c(0.05, 0.95)), 
                             chr=chr),
              double.bs=list(loci=these.loci,
                             pos=pos,
                             full.results=second.all.scans,
                             peak.loci=second.peak.loci.matrix,
                             peak.pos=second.peak.pos.list,
                             ci=quantile(second.peak.pos.list[[scale]], probs=c(0.05, 0.95)), 
                             chr=chr)))
}

################### Subsample CI
get.reduced.threshold.from.F <- function(original.threshold, N, k, prop){ # From Valdar et al. 2009 - based on the F test
  n <- floor(N*prop)
  return(-log10(pf(((n - k)/(N - k))*qf(p=10^-original.threshold, df1=k, df2=N, lower.tail=FALSE), df=k, df2=n, lower.tail=FALSE)))
}
subsample.ci <- function(formula, data, genomecache, K, model,
                         num.samples, num.av.over, num.qtl, sub.sample.prop=0.63, coverage.prob=c(0.95, 0.8), set.threshold=0,
                         num.imp, chr,
                         use.par="h2", fix.par=NULL, use.scan.fix.par=TRUE,
                         do.augment=FALSE, seed=1, rearrange=TRUE){
  
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/ci_support.R", local=TRUE)
  #source("/nas02/home/g/k/gkeele/SolbergHS/rats989/lmmbygls/GLS_03072016/support_functions.R", local=TRUE)
  
  full.to.dosages <- straineff.mapping.matrix()
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  
  cat("Beginning MI scans\n")
  
  loci <- h$getLoci()
  these.loci <- loci[h$getChromOfLocus(loci) == chr]
  
  ## Setting up data structures to record scan statistics
  peak.loci.matrix <- matrix(NA, nrow=num.samples, ncol=num.qtl)
  all.scans <- array(NA, dim=c(num.samples, num.qtl, length(these.loci)))
  
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula <- formula(paste0("y", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
  set.seed(seed)
  
  reduced.threshold <- get.reduced.threshold.from.F(original.threshold=set.threshold, N=nrow(data), k=7, prop=sub.sample.prop)
  
  sample.list <- list()
  count.successful.scans <- 0
  count.scans <- 0
  while(count.successful.scans < num.samples){
    this.subsample <- sort(sample.int(nrow(data), size=floor(sub.sample.prop*nrow(data)), replace=FALSE))
    pass <- FALSE # Set to FALSE, switch only gets flipped if QTL peaks pass threshold
    ## temporary storage
    subsample.scans <- matrix(NA, nrow=num.qtl, length(these.loci))
    peak.loci.vector <- rep(NA, num.qtl)
    for(j in 1:num.qtl){
      if(j == 1){
        this.data <- data[this.subsample,]
        this.K <- K[this.subsample, this.subsample]
        
        this.scan <- scan.h2lmm(genomecache=genomecache,
                                formula=this.formula,
                                data=this.data,
                                K=this.K,
                                model=model,
                                use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, brute=FALSE,
                                use.fix.par=use.scan.fix.par, chr=chr)
        pass <- ifelse(max(-log10(this.scan$p.value)) > reduced.threshold, TRUE, FALSE)
      }
      if(pass & j > 1){
        ## Probabilities of locus to regress out
        P <- h$getLocusMatrix(locus=peak.locus, model="full", subjects=this.data$SUBJECT.NAME)
        
        resid.vector <- regress.out.average.qtl()
        newer.data <- cbind(resid.vector, this.data)
        names(newer.data)[1] <- "resid"
        
        this.new.formula <- formula(paste0("resid", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
        this.scan <- scan.h2lmm(genomecache=genomecache,
                                formula=this.new.formula,
                                data=newer.data,
                                K=this.K,
                                model=model,
                                use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, brute=FALSE,
                                use.fix.par=use.scan.fix.par, chr=chr)
        pass <- ifelse(max(-log10(this.scan$p.value)) > reduced.threshold, TRUE, FALSE)
      }
      # Storing peak information after passes
      if(pass){
        peak.index <- which.max(-log10(this.scan$p.value))
        subsample.scans[j,] <- this.scan$p.value
        peak.locus <- this.scan$loci[peak.index]
        peak.loci.vector[j] <- peak.locus
        
        if(j == num.qtl){
          count.successful.scans <- count.successful.scans + 1
          
          all.scans[count.successful.scans,,] <- subsample.scans
          peak.loci.matrix[count.successful.scans,] <- peak.loci.vector
          sample.list[[count.successful.scans]] <- this.subsample
          
          cat("Finished scan of", count.successful.scans, "subsamples out of", num.samples, "\n")
        }
      }
      if(j == num.qtl){
        count.scans <- count.scans + 1
      }
      if(!pass & j == num.qtl){
        cat("Scan of", count.successful.scans + 1, "subsample out of", num.samples, "failed - ", count.scans, "samples attempted in total", "\n")
      }
    }
  }
  dimnames(all.scans)[[3]] <- loci[h$getChromOfLocus(loci) == chr]
  
  peak.cm.matrix <- t(apply(peak.loci.matrix, 1, function(x) h$getLocusStart(loci=x, scale="cM")))
  peak.mb.matrix <- t(apply(peak.loci.matrix, 1, function(x) h$getLocusStart(loci=x, scale="Mb")))
  
  # Reorder QTL
  if(rearrange & j > 1){
    order.these.positions <- peak.cm.matrix
    
    order.index <- apply(order.these.positions, 1, function(x) order(x))
    peak.cm.matrix <- t(sapply(1:num.samples, function(x) peak.cm.matrix[x, order.index[,x]]))
    peak.mb.matrix <- t(sapply(1:num.samples, function(x) peak.mb.matrix[x, order.index[,x]]))
    peak.loci.matrix <- t(sapply(1:num.samples, function(x) peak.loci.matrix[x, order.index[,x]]))
    for(i in 1:num.samples){
      all.scans[i,,] <- all.scans[i,order.index[,i],]
    }
  }
  peak.pos.list <- list(cM=peak.cm.matrix, Mb=peak.mb.matrix)
  pos <- list(cM=h$getLocusStart(these.loci, scale="cM"), Mb=h$getLocusStart(these.loci, scale="Mb"))
  ci <- list()
  for(i in 1:length(coverage.prob)){
    alpha <- 1 - coverage.prob[i]
    ci[[paste0("coverage", coverage.prob[i])]] <- list(cM=t(sapply(1:num.qtl, function(x) quantile(peak.pos.list[["cM"]][,x], probs=c(alpha/2, 1 - alpha/2), na.rm=TRUE))),
                                                       Mb=t(sapply(1:num.qtl, function(x) quantile(peak.pos.list[["Mb"]][,x], probs=c(alpha/2, 1 - alpha/2), na.rm=TRUE))))
  }
  return(list(loci=these.loci,
              pos=pos,
              peak.loci=peak.loci.matrix, 
              peak.pos=peak.pos.list,
              full.results=all.scans, 
              ci=ci, 
              chr=chr))
}

nonparametric.bootstrap <- function(formula, data, genomecache, K,
                                    num.bs.samples, num.imp, chr,
                                    model, use.par="h2", use.scan.fix.par=TRUE,
                                    do.augment=FALSE, seed=1, scale="cM"){
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  #source("/nas02/home/g/k/gkeele/SolbergHS/rats989/lmmbygls/GLS_03072016/support_functions.R", local=TRUE)
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  
  cat("Beginning MI scans\n")
  peak.loci.vec <- rep(NA, num.bs.samples)
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula <- formula(paste0("y", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
  set.seed(seed)
  
  weight.matrix <- matrix(NA, nrow=nrow(data), ncol=num.bs.samples)
  for(i in 1:num.bs.samples){
    x <- rep(0, nrow(data))
    names(x) <- 1:nrow(data)
    sample <- table(sort(sample.int(nrow(data), replace=TRUE)))
    x[names(sample)] <- 1/as.vector(sample)
    weight.matrix[,i] <- x
  }
  for(i in 1:num.bs.samples){
    nonzero <- which(weight.matrix[,i] != 0)
    these.weights <- weight.matrix[nonzero,i]
    this.data <- data[nonzero,]
    this.K <- K[nonzero, nonzero]
    this.bs.scan <- scan.h2lmm(genomecache=genomecache,
                               formula=this.formula,
                               data=this.data,
                               K=this.K,
                               model=model,
                               weights=these.weights,
                               use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, use.chol=TRUE, brute=FALSE,
                               use.fix.par=use.scan.fix.par, chr=chr)
    #browser()
    peak.index <- which.max(-log10(this.bs.scan$p.value))
    peak.loci.vec[i] <- this.bs.scan$loci[peak.index]
    
    cat("Finished scan of", i, "bootstrap sample out of", num.bs.samples, "\n")
  }
  peak.pos.vec <- h$getLocusStart(loci=peak.loci.vec, scale=scale)
  return(list(loci=peak.loci.vec, pos=peak.pos.vec, scale=scale, ci=quantile(peak.pos.vec, probs=c(0.05, 0.95)), chr=chr))
}

case.bootstrap <- function(formula, data, genomecache, K,
                           num.bs.samples, num.imp, chr,
                           model, use.par="h2", use.scan.fix.par=TRUE,
                           do.augment=FALSE, seed=1, scale="cM"){
  
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  #source("/nas02/home/g/k/gkeele/SolbergHS/rats989/lmmbygls/GLS_03072016/support_functions.R", local=TRUE)
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  
  cat("Beginning MI scans\n")
  peak.loci.vec <- rep(NA, num.bs.samples)
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula <- formula(paste0("y", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
  set.seed(seed)
  
  sample.list <- list()
  for(i in 1:num.bs.samples){
    sample.list[[i]] <- table(sort(sample.int(nrow(data), replace=TRUE)))
  }
  for(i in 1:num.bs.samples){
    resample.index <- rep(as.numeric(names(sample.list[[i]])), sample.list[[i]])
    this.data <- data[resample.index,]
    this.K <- K[resample.index, resample.index]
    this.bs.scan <- scan.h2lmm(genomecache=genomecache,
                               formula=this.formula,
                               data=this.data,
                               K=this.K,
                               model=model,
                               use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, use.chol=TRUE, brute=FALSE,
                               use.fix.par=use.scan.fix.par, chr=chr)
    peak.index <- which.max(-log10(this.bs.scan$p.value))
    peak.loci.vec[i] <- this.bs.scan$loci[peak.index]
    
    cat("Finished scan of", i, "bootstrap sample out of", num.bs.samples, "\n")
  }
  peak.pos.vec <- h$getLocusStart(loci=peak.loci.vec, scale=scale)
  return(list(loci=peak.loci.vec, pos=peak.pos.vec, scale=scale, ci=quantile(peak.pos.vec, probs=c(0.05, 0.95)), chr=chr))
}

bayesian.bootstrap <- function(formula, data, genomecache, K,
                               num.bs.samples, num.imp, chr,
                               model, use.par="h2", use.scan.fix.par=TRUE,
                               seed=1, scale="cM"){
  
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  #source("/nas02/home/g/k/gkeele/SolbergHS/rats989/lmmbygls/GLS_03072016/support_functions.R", local=TRUE)
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  
  set.seed(seed)
  
  cat("Generated", num.bs.samples, "bootstrap weights\n")
  cat("Beginning MI scans\n")
  
  peak.loci.vec <- rep(NA, num.bs.samples)
  
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula <- formula(paste0("y", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
  n <- nrow(data)
  
  weight.matrix <- matrix(NA, nrow=nrow(data), ncol=num.bs.samples)
  
  for(i in 1:num.bs.samples){
    random.num <- sort(runif(n-1, min = 0, max = 1))
    
    random.num[n] <- 1
    weights <- rep(NA, n)
    weights[1] <- random.num[1]
    for(j in 2:n){
      weights[j] <- random.num[j]-random.num[j-1]
    }
    weight.matrix[,i] <- 1/(weights*n)
  }
  for(i in 1:num.bs.samples){
    these.weights <- weight.matrix[,i]
    
    this.bs.scan <- scan.h2lmm(genomecache=genomecache,
                               formula=this.formula,
                               data=data,
                               K=K,
                               model=model, weights=these.weights,
                               use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, use.chol=TRUE, brute=FALSE,
                               use.fix.par=use.scan.fix.par, chr=chr)
    peak.index <- which.max(-log10(this.bs.scan$p.value))
    peak.loci.vec[i] <- this.bs.scan$loci[peak.index]
    
    cat("Finished scan of", i, "bootstrap sample out of", num.bs.samples, "\n")
  }
  peak.pos.vec <- h$getLocusStart(loci=peak.loci.vec, scale=scale)
  return(list(loci=peak.loci.vec, pos=peak.pos.vec, scale=scale, ci=quantile(peak.pos.vec, probs=c(0.05, 0.95)), chr=chr))
}

calc.LOD.drop <- function(scan.object, scale,
                          use.lod=TRUE, unit.drop=2){
  
  if(use.lod){ 
    outcome <- scan.object$LOD
  }
  if(!use.lod){ 
    outcome <- -log10(scan.object$p.value) 
  }
  pos <- unlist(scan.object$pos[scale])
  
  order.i <- order(pos)
  pos <- pos[order.i]
  outcome <- outcome[order.i]
  
  max.peak <- max(outcome)
  max.peak.index <- which.max(outcome)
  chr.of.locus <- scan.object$chr[max.peak.index]
  ub <- lb <- NULL
  move.up <- max.peak.index + 1
  move.down <- max.peak.index - 1
  while(is.null(ub)){
    if(outcome[move.up] < max.peak - unit.drop) {
      ub <- move.up
    }
    move.up <- move.up + 1
  }
  while(is.null(lb)){
    if(outcome[move.down] < max.peak - unit.drop) {
      lb <- move.down
    }
    move.down <- move.down - 1
  }
  lower.loci <- scan.object$loci[lb]
  upper.loci <- scan.object$loci[ub]
  return(list(scale=scale, loci.boundaries=c(lower.loci, upper.loci), ci=c(pos[lb], pos[ub]), chr=chr.of.locus))
}

nonparametric.bootstrap.deblup <- function(formula, data, genomecache, K,
                                           num.bs.samples, num.imp, chr,
                                           model, use.par="h2",
                                           do.augment=FALSE, seed=1, scale="cM"){
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/scan.h2lmm.deBLUP.R", local=TRUE)
  #source("/nas02/home/g/k/gkeele/SolbergHS/rats989/lmmbygls/GLS_03072016/support_functions.R", local=TRUE)
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  
  cat("Beginning MI scans\n")
  peak.loci.vec <- rep(NA, num.bs.samples)
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula <- formula(paste0("y", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
  set.seed(seed)
  
  weight.matrix <- matrix(NA, nrow=nrow(data), ncol=num.bs.samples)
  for(i in 1:num.bs.samples){
    x <- rep(0, nrow(data))
    names(x) <- 1:nrow(data)
    sample <- table(sort(sample.int(nrow(data), replace=TRUE)))
    x[names(sample)] <- 1/as.vector(sample)
    weight.matrix[,i] <- x
  }
  for(i in 1:num.bs.samples){
    nonzero <- which(weight.matrix[,i] != 0)
    these.weights <- weight.matrix[nonzero,i]
    this.data <- data[nonzero,]
    this.K <- K[nonzero, nonzero]
    this.bs.scan <- scan.h2lmm.deBLUP(genomecache=genomecache,
                                      formula=this.formula,
                                      data=this.data,
                                      K=this.K,
                                      model=model,
                                      weights=these.weights,
                                      use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, use.chol=TRUE, brute=FALSE,
                                      chr=chr)
    peak.index <- which.max(-log10(this.bs.scan$p.value))
    peak.loci.vec[i] <- this.bs.scan$loci[peak.index]
    
    cat("Finished scan of", i, "bootstrap sample out of", num.bs.samples, "\n")
  }
  peak.pos.vec <- h$getLocusStart(loci=peak.loci.vec, scale=scale)
  return(list(loci=peak.loci.vec, pos=peak.pos.vec, scale=scale, ci=quantile(peak.pos.vec, probs=c(0.05, 0.95)), chr=chr))
}

bayesian.bootstrap.deblup <- function(formula, data, genomecache, K,
                                      num.bs.samples, num.imp, chr,
                                      model, use.par="h2", use.scan.fix.par=TRUE,
                                      seed=1, scale="cM"){
  
  scan.environment <- environment()
  for(object in ls(scan.environment)){
    assign(object, get(object, scan.environment))
  }
  source("~/Documents/SolbergHS/GLS/lmmbygls/scripts_to_source.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/support_functions.R", local=TRUE)
  source("~/Documents/SolbergHS/GLS/lmmbygls/scan.h2lmm.deBLUP.R", local=TRUE)
  #source("/nas02/home/g/k/gkeele/SolbergHS/rats989/lmmbygls/GLS_03072016/support_functions.R", local=TRUE)
  
  h <- DiploprobReader$new(genomecache)
  cache.subjects <- h$getSubjects()
  
  data.and.K <- make.processed.data()
  data <- data.and.K$data
  K <- data.and.K$K
  
  set.seed(seed)
  
  cat("Generated", num.bs.samples, "bootstrap weights\n")
  cat("Beginning MI scans\n")
  
  peak.loci.vec <- rep(NA, num.bs.samples)
  
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula <- formula(paste0("y", " ~ ", unlist(strsplit(this.formula.string, split="~"))[-1]))
  n <- nrow(data)
  
  weight.matrix <- matrix(NA, nrow=nrow(data), ncol=num.bs.samples)
  
  for(i in 1:num.bs.samples){
    random.num <- sort(runif(n-1, min = 0, max = 1))
    
    random.num[n] <- 1
    weights <- rep(NA, n)
    weights[1] <- random.num[1]
    for(j in 2:n){
      weights[j] <- random.num[j]-random.num[j-1]
    }
    weight.matrix[,i] <- 1/(weights*n)
  }
  for(i in 1:num.bs.samples){
    these.weights <- weight.matrix[,i]
    this.bs.scan <- scan.h2lmm.deBLUP(genomecache=genomecache,
                                      formula=this.formula,
                                      data=data,
                                      K=K,
                                      model=model,
                                      weights=these.weights,
                                      use.par="h2", use.multi.impute=TRUE, num.imp=num.imp, use.chol=TRUE, brute=FALSE,
                                      chr=chr)
    peak.index <- which.max(-log10(this.bs.scan$p.value))
    peak.loci.vec[i] <- this.bs.scan$loci[peak.index]
    
    cat("Finished scan of", i, "bootstrap sample out of", num.bs.samples, "\n")
  }
  peak.pos.vec <- h$getLocusStart(loci=peak.loci.vec, scale=scale)
  return(list(loci=peak.loci.vec, pos=peak.pos.vec, scale=scale, ci=quantile(peak.pos.vec, probs=c(0.05, 0.95)), chr=chr))
}

