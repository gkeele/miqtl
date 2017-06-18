#' @export
null.pvalue.ci.plot <- function(null.scans, conf.level=0.95, scale="Mb",
                                main.label=NULL,
                                bs.max=NULL,
                                SIM.object=NULL,
                                dist.object=NULL,
                                cache.title=NULL,
                                formula.title=NULL,
                                model.title=NULL,
                                y.max.manual=NULL, y.min.manual=NULL,
                                give.output=FALSE,
                                is.sim=FALSE)
{
  chr <- null.scans$full.results$chr
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  
  pos <-null.scans$full.results$pos[[scale]]
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x))
  max.pos <- tapply(pos, pre.chr, function(x) max(x))
  chr.types <- levels(pre.chr)
  
  # Calculating confidence intervals
  these.cis <- apply(null.scans$full.results$p.values, 2, function(x) EnvStats::eexp(-log(x), ci=TRUE, conf.level=conf.level)$interval$limits[c("LCL", "UCL")])
  these.cis <- these.cis[,order.i]
  
  # Finding extremes of y for plot window
  y.max <- max(-log(these.cis[1,]), y.max.manual, 0.5)
  y.min <- min(-log(these.cis[2,]), y.min.manual, -0.5)
  
  # Adding color scale that matches siginficance
  if(is.null(SIM.object) & is.null(dist.object)){
    map.to.color <- as.vector(apply(abs(these.cis - 1), 2, function(x) min(x)))
    above.one <- (these.cis[1,] > 1) + (these.cis[2,] > 1)
    below.one <- (these.cis[1,] < 1) + (these.cis[2,] < 1)
    map.to.color[which(above.one == 2)] <- 0
    map.to.color[which(below.one == 2)] <- 0
    if(any(map.to.color != 0)){
      map.to.color <- floor((map.to.color/max(map.to.color))*1000)
    }
  }
  
  # Adding color based on information content
  if(!is.null(SIM.object)){
    map.to.color <- colMeans(SIM.object$SIM.matrix)
    map.to.color <- map.to.color[colnames(these.cis)]*1000
  }
  
  # Adding color based on imbalance of alleles
  if(!is.null(dist.object)){
    map.to.color <- dist.object
    map.to.color <- map.to.color[colnames(these.cis)]*1000
  } 

  gray.spectrum <- gray.colors(n=1000, start=0, end=0.8)
  lambda.code <- "\u03BB"
  
  null.formula.string <- ifelse(is.formula(null.scans$formula), Reduce(paste, deparse(null.scans$formula)), null.scans$formula)
  shift.left <- min(pos[pre.chr==chr.types[1]])
  if(!is.sim){
    this.title <- c(paste("%95 CI of -Log rate parameter", lambda.code, "from", nrow(null.scans$full.results$p.values), "null sims", sep=" "), 
                    cache.title,
                    paste0(null.formula.string, " + locus (", null.scans$model, ")")) 
  }
  if(is.sim){
    this.title <- c(paste("%95 CI of -Log rate parameter", lambda.code, "from", nrow(null.scans$full.results$p.values), "null sims", sep=" "), 
                    cache.title,
                    paste0("Parameters: n=", null.scans$num.sim, ", num.founders=", null.scans$num.founders))
  }

  plot(rep(pos[1], 2), -c(log(these.cis[1,1]), log(these.cis[2,1])), 
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)), type="l",
       ylim=c(y.min, y.max), 
       xaxt="n", xlab="", ylab=paste(conf.level, "CIs of -log", lambda.code), main=this.title,
       frame.plot=F, lwd=3, col=gray.spectrum[map.to.color[1]+1], las=1)
  if((these.cis[2,1] > 1 & these.cis[1,1] > 1) | (these.cis[2,1] < 1 & these.cis[1,1] < 1)){
    points(pos[1], y.max, pch=8, col="red")
  }

  this.pos <- pos[pre.chr == chr.types[1]]
  this.cis <- these.cis[, pre.chr == chr.types[1]]
  this.diff <- map.to.color[pre.chr == chr.types[1]]
  for(i in 2:length(this.pos)){
    #browser()
    lines(rep(this.pos[i], 2), -c(log(this.cis[1,i]), log(this.cis[2,i])), col=gray.spectrum[this.diff[i]+1])
    if((this.cis[2,i] > 1 & this.cis[1,i] > 1) | (this.cis[2,i] < 1 & this.cis[1,i] < 1)){
      points(this.pos[i], y.max, pch=8, col="red")
    }
  }
  
  car.blue <- "#7BAFD4"
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      this.cis <- these.cis[, pre.chr==chr.types[i]]
      this.diff <- map.to.color[pre.chr == chr.types[i]]
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos), min(this.pos):max(this.pos), max(this.pos)), 
                y=c(y.min, rep(y.max, length(min(this.pos):max(this.pos))), y.min), border=NA, col=car.blue)
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)
      for(j in 1:length(this.pos)){
        lines(rep(this.pos[j], 2), -c(log(this.cis[1,j]), log(this.cis[2,j])), col=gray.spectrum[this.diff[j]+1])
        #points(rep(this.pos[j], 2), c(this.cis[1,j], this.cis[2,j]), pch="-", col=gray.spectrum[this.diff[j]], cex=1)
        if((this.cis[2,j] > 1 & this.cis[1,j] > 1) | (this.cis[2,j] < 1 & this.cis[1,j] < 1)){
          points(this.pos[j], y.max, pch=8, col="red")
        }
      }
      shift <- shift + max.pos[i]
    }
  }
  abline(h=0, col="red", lty=2)
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  if(!has.X){
    axis.label <- chr.types
  }
  axis(side=1, tick=FALSE, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-3.5)
  if(give.output){
    return(these.cis)
  }
}


null.pvalue.plot <- function(par.bs.scans, conf.level=0.95,
                             main.label=NULL,
                             bs.max=NULL,
                             cache.title="DO-QTL",
                             formula.title=NULL,
                             model.title=NULL,
                             y.max.manual=NULL, y.min.manual=NULL, 
                             give.output=FALSE,
                             is.sim=FALSE)
{
  chr <- par.bs.scans$full.results$chr
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  
  pos <-par.bs.scans$full.results$pos
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x))
  max.pos <- tapply(pos, pre.chr, function(x) max(x))
  chr.types <- levels(pre.chr)
  
  num.samples <- nrow(par.bs.scans$full.results$p.values)
  
  
  # Confidence intervals
  these.cis <- apply(par.bs.scans$full.results$p.values, 2, function(x) EnvStats::eexp(-log(x), ci=TRUE, conf.level=conf.level)$interval$limits[c("LCL", "UCL")])
  these.cis <- these.cis[,order.i]
  
  midpoint.ci <- colMeans(these.cis)
  
  # Calculating p-values
  rate.mle <- 1/colMeans(-log(par.bs.scans$full.results$p.values))
  logLik1 <- sapply(1:length(rate.mle), function(x) num.samples*log(rate.mle[x]) - num.samples)
  logLik0 <- sapply(1:length(rate.mle), function(x) - num.samples/rate.mle[x])
  pval <- pchisq(q=-2*(logLik0 - logLik1), df=1, lower.tail=FALSE)
  pval[pval == 0] <- .Machine$double.xmin
  pval <- pval[order.i]
  logp <- -log10(pval)
  logp[midpoint.ci > 1] <- logp[midpoint.ci > 1] * -1
  
  # Finding extremes of y for plot window
  y.max <- ifelse(length(logp[logp > 0]) > 0, max(logp[logp > 0], y.max.manual), 1)
  y.min <- ifelse(length(logp[logp < 0]) > 0, min(logp[logp < 0], y.min.manual), -1)
  
  shift.left <- min(pos[chr==chr.types[1]])
  fpr <- mean(pval <= 0.05)
  lambda.code <- "\u03BB"
  if(!is.sim){
    this.title <- c(paste("Test of rate paramter", lambda.code, "for", nrow(par.bs.scans$full.results$p.values), "null sims"),
                    cache.title,
                    paste0("FPR=", round(fpr, 5))) 
  }
  if(is.sim){
    this.title <- c(paste("Test of rate paramter", lambda.code, "for", nrow(par.bs.scans$full.results$p.values), "null sims"),
                    cache.title,
                    paste0("Parameters: n=", par.bs.scans$num.sim, ", num.founders=", par.bs.scans$num.founders, ", FPR=", round(fpr, 5)))
  }
  plot(c(pos[1], pos[1]), c(0, logp[1]), 
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)), type="l",
       ylim=c(y.min, y.max), 
       xaxt="n", xlab="", ylab=expression('-log'[10]*'P'), main=this.title,
       frame.plot=F, lwd=3, col="gray", las=1)

  this.pos <- pos[pre.chr == chr.types[1]]
  this.logp <- logp[pre.chr == chr.types[1]]
  for(i in 2:length(this.pos)){
    lines(rep(this.pos[i], 2), c(0, this.logp[i]), col="black")
  }
  
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      this.logp <- logp[pre.chr==chr.types[i]]
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos), min(this.pos):max(this.pos), max(this.pos)), 
                y=c(y.min, rep(y.max, length(min(this.pos):max(this.pos))), y.min), border=NA, col="gray")
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)
      for(j in 1:length(this.pos)){
        lines(rep(this.pos[j], 2), c(0, this.logp[j]), col="black")
      }
      shift <- shift + max.pos[i]
    }
  }
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  if(!has.X){
    axis.label <- chr.types
  }
  axis(side=1, tick=F, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-3.5)
  if(give.output){
    return(pval)
  }
}


