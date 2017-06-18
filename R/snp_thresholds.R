#' @export
snp.null.par.bs.threshold.scan <- function(scan.object, 
                                           data,
                                           allele.dir, genomecache,
                                           model=c("additive", "full"), chr="all", use.REML=TRUE,
                                           use.par="h2", brute=TRUE, use.fix.par=FALSE, seed=1,
                                           use.chol=FALSE,
                                           just.these.loci=NULL,
                                           print.locus.fit=FALSE,
                                           map.file=NULL,
                                           num.bs.scans=100, keep.full.scans=TRUE,
                                           ...){
  h <- DiploprobReader$new(genomecache)
  
  loci <- scan.object$loci
  loci.chr <- scan.object$chr
  
  full.results <- these.chr <- these.pos <- NULL
  if(keep.full.scans){
    full.results <- matrix(NA, nrow=num.bs.scans, ncol=length(loci))
    these.pos <- scan.object$pos
  }
  max.results <- rep(NA, num.bs.scans)
  
  if(use.REML){
    null.fit <- scan.object$fit0.REML
  }
  if(!use.REML){
    null.fit <- scan.object$fit0
  }
  Xb <- null.fit$x %*% null.fit$coefficients
  
  set.seed(seed)
  new.y.matrix <- matrix(NA, nrow=length(Xb), ncol=num.bs.scans)
  for(i in 1:num.bs.scans){
    u <- c(mnormt::rmnorm(1, mean=rep(0, nrow(null.fit$x)), varcov=null.fit$K*null.fit$tau2.mle))
    e <- rnorm(n=nrow(null.fit$x), mean=0, sd=sqrt(null.fit$sigma2.mle))
    new.y.matrix[,i] <- Xb + u + e
  }
  for(i in 1:num.bs.scans){
    new.y <- data.frame(new.y=new.y.matrix[,i], SUBJECT.NAME=colnames(null.fit$K))
    par.bs.data <- merge(x=new.y, y=data, by="SUBJECT.NAME", all.x=TRUE)
    
    original.formula <- paste0(Reduce(paste, deparse(formula(null.fit))))
    par.bs.formula <- formula(paste0("new.y ~ ", unlist(strsplit(original.formula, split="~"))[-1]))
    cat(paste("Conducting null bootstrap scan:", i, "\n"))
    par.bs.scan <- imputed.snp.scan.h2lmm(data=par.bs.data, formula=par.bs.formula, K=null.fit$K,
                                          allele.dir=allele.dir, genomecache=genomecache, 
                                          model=model,
                                          use.par=use.par, chr=chr, brute=brute, use.fix.par=use.fix.par, seed=seed, 
                                          use.chol=use.chol,
                                          just.these.loci=NULL, print.locus.fit=print.locus.fit,
                                          X.list=scan.object$X.list, return.X.list=FALSE)
    if(keep.full.scans){
      full.results[i,] <- par.bs.scan$p.value
    }
    max.results[i] <- min(par.bs.scan$p.value)
  }
  return(list(full.results=list(p.values=full.results, chr=these.chr, pos=these.pos), max.statistics=max.results))
}

#' @export
snp.par.perm.threshold.scan <- function(scan.object, 
                                        data,
                                        allele.dir, genomecache,
                                        model=c("additive", "full"), chr="all", use.REML=TRUE,
                                        use.par="h2", brute=TRUE, use.fix.par=FALSE, seed=1,
                                        use.chol=FALSE,
                                        just.these.loci=NULL,
                                        print.locus.fit=FALSE,
                                        map.file=NULL,
                                        num.perm.scans=100, keep.full.scans=TRUE,
                                        ...){
  h <- DiploprobReader$new(genomecache)
  
  loci <- scan.object$loci
  loci.chr <- scan.object$chr
  
  full.results <- these.chr <- these.pos <- NULL
  if(keep.full.scans){
    full.results <- matrix(NA, nrow=num.perm.scans, ncol=length(loci))
    these.pos <- scan.object$pos
  }
  max.results <- rep(NA, num.perm.scans)
  
  if(use.REML){
    null.fit <- scan.object$fit0.REML
  }
  if(!use.REML){
    null.fit <- scan.object$fit0
  }
  Xb <- null.fit$x %*% null.fit$coefficients
  
  set.seed(seed)
  perm.y.matrix <- matrix(NA, nrow=length(Xb), ncol=num.perm.scans)
  for(i in 1:num.perm.scans){
    u <- c(mnormt::rmnorm(1, mean=rep(0, nrow(null.fit$x)), varcov=null.fit$K*null.fit$tau2.mle))
    e <- rnorm(n=nrow(null.fit$x), mean=0, sd=sqrt(null.fit$sigma2.mle))
    perm.y.ranks <- order(Xb + u + e)
    perm.y.matrix[,i] <- null.fit$y[perm.y.ranks]
  }
  for(i in 1:num.perm.scans){
    perm.y <- data.frame(perm.y=perm.y.matrix[,i], SUBJECT.NAME=colnames(null.fit$K))
    perm.data <- merge(x=perm.y, y=data, by="SUBJECT.NAME", all.x=TRUE)
    
    original.formula <- paste0(Reduce(paste, deparse(formula(null.fit))))
    perm.formula <- formula(paste0("perm.y ~ ", unlist(strsplit(original.formula, split="~"))[-1]))
    cat(paste("Conducting permutation scan:", i, "\n"))
    perm.scan <- imputed.snp.scan.h2lmm(data=perm.data, formula=perm.formula, K=null.fit$K,
                                        allele.dir=allele.dir, genomecache=genomecache, 
                                        model=model,
                                        use.par=use.par, chr=chr, brute=brute, use.fix.par=use.fix.par, seed=seed, 
                                        use.chol=use.chol,
                                        just.these.loci=NULL, print.locus.fit=print.locus.fit,
                                        X.list=scan.object$X.list, return.X.list=FALSE)
    if(keep.full.scans){
      full.results[i,] <- perm.scan$p.value
    }
    max.results[i] <- min(perm.scan$p.value)
  }
  return(list(full.results=list(p.values=full.results, chr=these.chr, pos=these.pos), max.statistics=max.results))
}

#' @export
snp.perm.threshold.scan <- function(scan.object, 
                                    data,
                                    allele.dir, genomecache,
                                    model=c("additive", "full"), chr="all",
                                    use.par="h2", brute=TRUE, use.fix.par=FALSE, seed=1,
                                    use.chol=FALSE,
                                    just.these.loci=NULL,
                                    print.locus.fit=FALSE,
                                    map.file=NULL,
                                    num.perm.scans=100, keep.full.scans=TRUE,
                                    ...){
  h <- DiploprobReader$new(genomecache)
  
  loci <- scan.object$loci
  loci.chr <- scan.object$chr
  
  full.results <- these.chr <- these.pos <- NULL
  if(keep.full.scans){
    full.results <- matrix(NA, nrow=num.perm.scans, ncol=length(loci))
    these.pos <- scan.object$pos
  }
  max.results <- rep(NA, num.perm.scans)
  
  set.seed(seed)
  perm.y.matrix <- matrix(NA, nrow=length(scan.object$fit0$y), ncol=num.perm.scans)
  for(i in 1:num.perm.scans){
    perm.y.matrix[,i] <- sample(1:nrow(perm.y.matrix))
  }
  for(i in 1:num.perm.scans){
    perm.y <- data.frame(perm.y=perm.y.matrix[,i], SUBJECT.NAME=colnames(scan.object$fit0$K))
    perm.data <- merge(x=perm.y, y=data, by="SUBJECT.NAME", all.x=TRUE)
    
    original.formula <- paste0(Reduce(paste, deparse(formula(scan.object$formula))))
    perm.formula <- formula(paste0("perm.y ~ ", unlist(strsplit(original.formula, split="~"))[-1]))
    cat(paste("Conducting permutation scan:", i, "\n"))
    perm.scan <- imputed.snp.scan.h2lmm(data=perm.data, formula=perm.formula, K=scan.object$fit0$K,
                                        allele.dir=allele.dir, genomecache=genomecache, 
                                        model=model,
                                        use.par=use.par, chr=chr, brute=brute, use.fix.par=use.fix.par, seed=seed, 
                                        use.chol=use.chol,
                                        just.these.loci=NULL, print.locus.fit=print.locus.fit,
                                        X.list=scan.object$X.list, return.X.list=FALSE)
    if(keep.full.scans){
      full.results[i,] <- perm.scan$p.value
    }
    max.results[i] <- min(perm.scan$p.value)
  }
  return(list(full.results=list(p.values=full.results, chr=these.chr, pos=these.pos), max.statistics=max.results))
}