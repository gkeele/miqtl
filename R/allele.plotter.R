#' Plot the regression coefficients or BLUPs as allele effects for a genome scan (whole genomes or subset of chromosomes)
#'
#' This function takes a genome scan object from scan.h2lmm() and plots out allele effect estimates
#' based on regression coefficients or BLUPs.
#'
#' @param scan.object A scan.h2lmm() object that had return.allele.effects=TRUE.
#' @param just.these.chr DEFAULT: NULL. Specifies a subset of the chromosomes to be plotted. NULL results in all chromosomes being plotted.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param main.colors DEFAULT: CC founder colors. The colors that correspond to the founder alleles.
#' @param use.legend DEFAULT: TRUE. Include a legend for the different associations. If TRUE, the labels are the names of the non.mi.scan.list object.
#' @param main DEFAULT: NULL. Adds a title above the model.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param my.lwd DEFAULT: rep(1.5, 8). Specifies the lwds for the alleles.
#' @param my.legend.pos DEFAULT: "topright". Specifies where to put the legend, if specified in use.legend.
#' @param transparency DEFAULT: 0.6. The transparency level for the confidence bands. 0 is invisible, 1 is solid.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. 
#' @param y.min.manual DEFAULT: NULL. Manually adds a min y-value.
#' @param no.title DEFAULT: FALSE. If TRUE, produces a plot with no title, not even the automatic
#' title components like formula and sample size.
#' @param override.title DEFAULT: NULL. If non-null, overrides automatic title components like
#' formula and sample size.
#' @param add.chr.to.label DEFAULT: FALSE. If TRUE, "chr" precedes integers on the x-axis.
#' @param alternative.labels DEFAULT: NULL. If non-null, specifies alternative labels for the
#' founder alleles.
#' @param axis.cram DEFAULT: TRUE. This makes the plot much more likely to include all chromosome labels. With small plots, this could
#' lead to overlapping labels.
#' @param include.x.axis.line DEFAULT: TRUE. IF TRUE, this option adds an x-axis line with ticks between chromosomes.
#' @export
#' @examples allele.plotter.whole()
allele.plotter.whole <- function(scan.object, just.these.chr=NULL,
                                 scale="Mb", imp.confint.alpha=NULL,
                                 main.colors=c(rgb(240, 240, 0, maxColorValue=255), # yellow
                                               rgb(128, 128, 128, maxColorValue=255), # gray
                                               rgb(240, 128, 128, maxColorValue=255), # pink/salmon
                                               rgb(16, 16, 240, maxColorValue=255), # dark blue
                                               rgb(0, 160, 240, maxColorValue=255), # light blue
                                               rgb(0, 160, 0, maxColorValue=255), # green
                                               rgb(240, 0, 0, maxColorValue=255), # red
                                               rgb(144, 0, 224, maxColorValue=255)), # purple,
                                 use.legend=TRUE, main="", my.bty="n", my.lwd=1.25,
                                 set.plot.limit=c(-10, 10), # Null places no limit on y-axis
                                 my.legend.cex=0.6, my.legend.pos="topright", transparency=0.6,
                                 y.max.manual=NULL, y.min.manual=NULL, my.y.line=2, my.y.axis.cex=1, my.y.lab.cex=0.7,
                                 my.x.lab.cex=0.7,
                                 no.title=FALSE, override.title=NULL,
                                 add.chr.to.label=FALSE, alternative.labels=NULL,
                                 axis.cram=TRUE, include.x.axis.line=TRUE)
{
  allele.effects <- scan.object$allele.effects
  
  if(length(my.lwd) == 1){
    my.lwd <- rep(my.lwd, 8)
  }
  
  ## Check that scan object has allele effect estimate - stop if not
  if(is.null(allele.effects)){
    stop("No allele effects in scan object. Re-run scan.h2lmm with return.allele.effects=TRUE", call.=FALSE)
  }
  ## Check that confint is specified, imputations were used
  if(length(dim(allele.effects)) != 3 & !is.null(imp.confint.alpha)){
    cat("Multiple imputations must have been used to plot effect confidence intervals!\n",
        "Setting imp.confint.alpha to NULL\n")
    imp.confint.alpha <- NULL
  }
  
  chr <- scan.object$chr
  pos <- ifelse(rep(scale=="Mb", length(scan.object$loci)), scan.object$pos$Mb, scan.object$pos$cM)

  ## Grabbing the number of alleles
  if(length(dim(allele.effects)) == 3){ num.founders <- dim(allele.effects)[1] }
  else{ num.founders <- nrow(allele.effects) }
    
  if(!is.null(just.these.chr)){
    keep.chr <- chr %in% just.these.chr
    chr <- chr[keep.chr]
    if(length(dim(allele.effects)) == 3){
      allele.effects <- allele.effects[, keep.chr,]
    }
    else{
      allele.effects <- allele.effects[, keep.chr]
    }
    pos <- pos[keep.chr]
  }
  
  ## Removing NAs that start in pos
  pos.na <- is.na(pos)
  pos <- pos[!(pos.na)]
  chr <- chr[!(pos.na)]
  
  ## Ordering loci
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  if(length(dim(allele.effects)) == 3){
    allele.effects <- allele.effects[, order.i,]
  }
  else{
    allele.effects <- allele.effects[, order.i]
  }
  
  ## Calculating confidence intervals on the allele effect means
  allele.effect.confint <- NULL
  if(!is.null(imp.confint.alpha)){
    allele.effect.confint <- apply(allele.effects, c(1, 2), function(x) ci.mean(x, alpha=imp.confint.alpha))
  }
  ## Calculating the allele effect means given that multiple imputations were used
  if(length(dim(allele.effects)) == 3){
    allele.effects <- apply(allele.effects, c(1, 2), function(x) mean(x))
  }
  
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- max(as.numeric(unique(chr[chr != "X"]))) + 1
  }
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x, na.rm=TRUE))
  max.pos <- tapply(pos, pre.chr, function(x) max(x, na.rm=TRUE))
  chr.types <- levels(pre.chr)
  
  # Finding max and min of y for plot window
  y.max <- ifelse(max(allele.effects, na.rm=TRUE) > 0, 
                  ceiling(max(allele.effects, allele.effect.confint, na.rm=TRUE)), 
                  floor(max(allele.effects, allele.effect.confint, na.rm=TRUE)))
  y.min <- ifelse(min(allele.effects, na.rm=TRUE) > 0, 
                  ceiling(min(allele.effects, allele.effect.confint, na.rm=TRUE)), 
                  floor(min(allele.effects, allele.effect.confint, na.rm=TRUE)))
  
  if(!is.null(set.plot.limit)){
    y.max <- min(y.max, set.plot.limit[2])
    y.min <- max(y.min, set.plot.limit[1])
  }
  
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  if(!is.null(y.min.manual)){
    y.min <- y.min.manual
  }
  shift.left <- min(pos[chr==chr.types[1]], na.rm=TRUE)
  
  ### Fixef or ranef
  locus.effect.type <- ifelse(scan.object$locus.effect.type == "fixed", "fixef", "ranef")
  locus.term <- paste("locus", locus.effect.type, sep=".")
  
  ### Handling the annoying differences between lmer and lm objects
  if(is.null(scan.object$fit0)){
    this.title <- c(main, 
                    paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"),
                    paste("n =", round(scan.object$n, 2)))
  }
  else{
    if(class(scan.object$fit0) != "lmerMod"){
      this.title <- c(main, 
                      paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"),
                      paste("n =", round(ifelse(is.null(scan.object$fit0$weights), 
                                                length(scan.object$fit0$y),
                                                sum(scan.object$fit0$weights)), 2)))
    }
    else{
      this.title <- c(main, 
                      paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"),
                      paste("n =", round(sum(scan.object$fit0@resp$weights), 2)))
    }
  }
  if(no.title){ this.title <- NULL }
  if(!is.null(override.title)){ this.title <- override.title }
  
  x.max <- sum(max.pos)+(length(chr.types)-1)
  plot(0, pch="",
       xlim=c(0, x.max), 
       ylim=c(y.min, y.max), 
       xaxt="n", yaxt="n", xlab="", ylab="", main=this.title,
       frame.plot=FALSE, type="l", cex=0.5, lwd=1.5, col=main.colors[1])
  axis(side=2, at=sort(unique(c(0, y.min, y.min:y.max, y.max))), las=2, cex.axis=my.y.axis.cex)
  mtext(text="Additive allele effects", side=2, line=my.y.line, cex=my.y.lab.cex)
  label.spots <- max.pos[1]/2
  x.tick.spots <- c(0, max.pos[1])
  
  ## Confint for first chr
  if(!is.null(imp.confint.alpha)){
    first.chr.x <- pos[pre.chr==chr.types[1]]
    for(i in 1:num.founders){
      founder.first.chr.y.max <- allele.effect.confint[2, i,  pre.chr==chr.types[1]]
      founder.first.chr.y.min <- allele.effect.confint[1, i,  pre.chr==chr.types[1]]
      if(any(is.na(founder.first.chr.y.max))){
        polygon.with.nas(x=first.chr.x,
                         y.ci.max=founder.first.chr.y.max, y.ci.min=founder.first.chr.y.min,
                         col=scales::alpha(main.colors[i], transparency))
      }
      else{
        polygon(x=c(first.chr.x, rev(first.chr.x)),
                y=c(founder.first.chr.y.max, rev(founder.first.chr.y.min)),
                border=NA, col=scales::alpha(main.colors[i], transparency))
      }
    }
  }
  ## Means for first chr
  for(i in 1:num.founders){
    points(x=pos[pre.chr==chr.types[1]], y=allele.effects[i,  pre.chr==chr.types[1]],
           type="l", lwd=my.lwd[i], col=main.colors[i])
  }
  if(length(chr.types) > 1){
    shift <- max.pos[1]
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        polygon(x=c(shift, shift:max(this.pos, na.rm=TRUE), max(this.pos, na.rm=TRUE)), 
                y=c(y.min, rep(y.max, length(shift:max(this.pos, na.rm=TRUE))), y.min), border=NA, col="gray88")
      }
      
      ## Confint for later chr
      if(!is.null(imp.confint.alpha)){
        this.chr.x <- pos[pre.chr==chr.types[i]] + shift
        for(j in 1:num.founders){
          founder.this.chr.y.max <- allele.effect.confint[2, j,  pre.chr==chr.types[i]]
          founder.this.chr.y.min <- allele.effect.confint[1, j,  pre.chr==chr.types[i]]
          if(any(is.na(founder.this.chr.y.max))){
            polygon.with.nas(x=this.chr.x,
                             y.ci.max=founder.this.chr.y.max, y.ci.min=founder.this.chr.y.min,
                             col=scales::alpha(main.colors[j], transparency))
          }
          else{
            polygon(x=c(this.chr.x, rev(this.chr.x)),
                    y=c(founder.this.chr.y.max, rev(founder.this.chr.y.min)),
                    border=NA, col=scales::alpha(main.colors[j], transparency))
          }
        }
      }
      ## Means for later chr
      for(j in 1:num.founders){
        points(x=pos[pre.chr==chr.types[i]] + shift, y=allele.effects[j,  pre.chr==chr.types[i]],
               type="l", lwd=my.lwd[j], col=main.colors[j])
      }
      label.spots <- c(label.spots, shift + max.pos[i]/2)
      x.tick.spots <- c(x.tick.spots, max.pos[i] + shift)
      shift <- shift + max.pos[i]
    }
  }
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  else{
    axis.label <- chr.types
  }
  
  if(include.x.axis.line){
    axis(side=1, tick=TRUE, line=NA, at=x.tick.spots, 
         labels=NA, xpd=TRUE, )
  }
  
  if(add.chr.to.label){
    axis.label <- paste("Chr", axis.label)
  }
  else{
    axis.label <- c("Chr", axis.label)
    label.spots <- c(-0.04*x.max, label.spots)
  }
  
  if(axis.cram){
    odd.axis.label <- axis.label[(1:length(axis.label) %% 2) == 1]
    odd.label.spots <- label.spots[(1:length(label.spots) %% 2) == 1]
    
    even.axis.label <- axis.label[(1:length(axis.label) %% 2) == 0]
    even.label.spots <- label.spots[(1:length(label.spots) %% 2) == 0]
    
    axis(side=1, tick=FALSE, line=NA, at=odd.label.spots, labels=odd.axis.label, cex.axis=my.x.lab.cex, padj=-1.5, xpd=TRUE)
    axis(side=1, tick=FALSE, line=NA, at=even.label.spots, labels=even.axis.label, cex.axis=my.x.lab.cex, padj=-1.5, xpd=TRUE)
  }
  else{
    axis(side=1, tick=FALSE, line=NA, at=label.spots, labels=axis.label, cex.axis=my.x.lab.cex, padj=-1.5, xpd=TRUE)
  }
  if(use.legend){
    if(!is.null(alternative.labels)){ these.labels <- alternative.labels }
    else{ these.labels <- rownames(allele.effects) }
    legend(my.legend.pos, legend=these.labels, 
           fill=main.colors, border="black", 
           col=main.colors, bty=my.bty, cex=my.legend.cex, bg="white")
  }
}

#' @export
allele.plotter.region <- function(scan.object,
                                  scale="Mb", 
                                  chr, region.min=NULL, region.max=NULL,
                                  imp.confint.alpha=NULL,
                                  main.colors=c(rgb(240, 240, 0, maxColorValue=255), # yellow
                                               rgb(128, 128, 128, maxColorValue=255), # gray
                                               rgb(240, 128, 128, maxColorValue=255), # pink/salmon
                                               rgb(16, 16, 240, maxColorValue=255), # dark blue
                                               rgb(0, 160, 240, maxColorValue=255), # light blue
                                               rgb(0, 160, 0, maxColorValue=255), # green
                                               rgb(240, 0, 0, maxColorValue=255), # red
                                               rgb(144, 0, 224, maxColorValue=255)), # purple,
                                  use.legend=TRUE, main="", my.bty="n", my.lwd=rep(1.25, 8),
                                  set.plot.limit=c(-5, 5), # Null places no limit on y-axis
                                  my.legend.cex=0.7, my.legend.pos="topright", transparency=0.6,
                                  y.max.manual=NULL, y.min.manual=NULL, 
                                  my.y.line=2, my.y.axis.cex=1, my.y.lab.cex=0.5,
                                  my.title.line=0.5, my.title.cex=1,
                                  no.title=FALSE, override.title=NULL,
                                  alternative.labels=NULL,
                                  rug.pos=NULL, rug.col="gray50")
{
  allele.effects <- scan.object$allele.effects
  
  if(length(my.lwd) == 1){
    my.lwd <- rep(my.lwd, 8)
  }
  
  ## Check that scan object has allele effect estimate - stop if not
  if(is.null(allele.effects)){
    stop("No allele effects in scan object. Re-run scan.h2lmm with return.allele.effects=TRUE", call.=FALSE)
  }
  ## Check that confint is specified, imputations were used
  if(length(dim(allele.effects)) != 3 & !is.null(imp.confint.alpha)){
    cat("Multiple imputations must have been used to plot effect confidence intervals!\n",
        "Setting imp.confint.alpha to NULL\n")
    imp.confint.alpha <- NULL
  }
  
  pos <- ifelse(rep(scale=="Mb", length(scan.object$loci)), scan.object$pos$Mb, scan.object$pos$cM)
  
  ## Getting index to plot
  keep.chr <- scan.object$chr == chr
  keep.na <- !is.na(pos)
  keep <- (keep.chr + keep.na) == 2
  order.i <- order(pos[keep])
  pos <- pos[keep][order.i]
  
  ## Grabbing the number of alleles
  if(length(dim(allele.effects)) == 3){ 
    num.founders <- dim(allele.effects)[1]
    allele.effects <- allele.effects[,keep,][,order.i,]
  }
  else{ 
    num.founders <- nrow(allele.effects) 
    allele.effects <- allele.effects[,keep][,order.i]
  }
  
  ## Calculating confidence intervals on the allele effect means
  allele.effect.confint <- NULL
  if(!is.null(imp.confint.alpha)){
    allele.effect.confint <- apply(allele.effects, c(1, 2), function(x) ci.mean(x, alpha=imp.confint.alpha))
  }
  ## Calculating the allele effect means given that multiple imputations were used
  if(length(dim(allele.effects)) == 3){
    allele.effects <- apply(allele.effects, c(1, 2), function(x) mean(x))
  }
  
  min.pos <- min(pos)
  max.pos <- max(pos)
  if(!is.null(region.min)){ min.pos <- region.min }
  if(!is.null(region.max)){ max.pos <- region.max }

  # Finding max and min of y for plot window
  y.max <- ifelse(max(allele.effects, na.rm=TRUE) > 0, 
                  ceiling(max(allele.effects, allele.effect.confint, na.rm=TRUE)), 
                  floor(max(allele.effects, allele.effect.confint, na.rm=TRUE)))
  y.min <- ifelse(min(allele.effects, na.rm=TRUE) > 0, 
                  ceiling(min(allele.effects, allele.effect.confint, na.rm=TRUE)), 
                  floor(min(allele.effects, allele.effect.confint, na.rm=TRUE)))
  
  if(!is.null(set.plot.limit)){
    y.max <- min(y.max, set.plot.limit[2])
    y.min <- max(y.min, set.plot.limit[1])
  }
  
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  if(!is.null(y.min.manual)){
    y.min <- y.min.manual
  }

  ### Fixef or ranef
  locus.effect.type <- ifelse(scan.object$locus.effect.type == "fixed", "fixef", "ranef")
  locus.term <- paste("locus", locus.effect.type, sep=".")
  
  ### Handling the annoying differences between lmer and lm objects
  if(is.null(scan.object$fit0)){
    this.title <- c(main, 
                    paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"),
                    paste("n =", round(scan.object$n, 2)))
  }
  else{
    if(class(scan.object$fit0) != "lmerMod"){
      this.title <- c(main, 
                      paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"),
                      paste("n =", round(ifelse(is.null(scan.object$fit0$weights), 
                                                length(scan.object$fit0$y),
                                                sum(scan.object$fit0$weights)), 2)))
    }
    else{
      this.title <- c(main, 
                      paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"),
                      paste("n =", round(sum(scan.object$fit0@resp$weights), 2)))
    }
  }
  if(no.title){ this.title <- NULL }
  if(!is.null(override.title)){ this.title <- override.title }
  
  this.xlab <- paste0("Chr ", chr, " Position (", scale, ")")
  
  plot(0, pch="",
       xlim=c(min.pos, max.pos), 
       ylim=c(y.min, y.max), 
       xaxt="n", yaxt="n", xlab=this.xlab, ylab="", main="",
       frame.plot=FALSE, type="l", cex=0.5, lwd=1.5)
  title(main=this.title, line=my.title.line, cex.main=my.title.cex)
  label.spots <- max.pos[1]/2
  x.tick.spots <- c(0, max.pos[1])
  
  x.ticks <- seq(min.pos, max.pos, length.out=5)
  x.ticks <- round(x.ticks)
  axis(side=1, tick=TRUE, line=NA, at=x.ticks, xpd=TRUE)
  axis(side=2, at=sort(unique(c(0, y.min, y.min:y.max, y.max))), las=2, cex.axis=my.y.axis.cex)
  mtext(text="Additive allele effects", side=2, line=my.y.line, cex=my.y.lab.cex)
  
  ## Confint
  if(!is.null(imp.confint.alpha)){
    for(i in 1:num.founders){
      founder.y.max <- allele.effect.confint[2, i,]
      founder.y.min <- allele.effect.confint[1, i,]
      if(any(is.na(founder.y.max))){
        polygon.with.nas(x=pos,
                         y.ci.max=founder.y.max, y.ci.min=founder.y.min,
                         col=scales::alpha(main.colors[i], transparency))
      }
      else{
        polygon(x=c(pos, rev(pos)),
                y=c(founder.y.max, rev(founder.y.min)),
                border=NA, col=scales::alpha(main.colors[i], transparency))
      }
    }
  }
  ## Means
  for(i in 1:num.founders){
    points(x=pos, y=allele.effects[i,],
           type="l", lwd=my.lwd[i], col=main.colors[i])
  }
  
  if(use.legend){
    if(!is.null(alternative.labels)){ these.labels <- alternative.labels }
    else{ these.labels <- rownames(allele.effects) }
    legend(my.legend.pos, legend=these.labels, 
           fill=main.colors, border="black", 
           col=main.colors, bty=my.bty, bg="white", cex=my.legend.cex)
  }
  if(!is.null(rug.pos)){
    if(length(rug.col) == 1){ rug.col <- rep(rug.col, length(rug.pos))}
    for(i in 1:length(rug.pos)){
      axis(1, at=rug.pos[i], col.ticks=rug.col[i], label=FALSE, cex.axis=0.1, lwd.ticks=3, lend="butt", tck=-0.1)
    }
  }
}

polygon.with.nas <- function(x, y.ci.max, y.ci.min, col){
  enc <- rle(!is.na(y.ci.max))
  endIdxs <- cumsum(enc$lengths)
  for(i in 1:length(enc$lengths)){
    if(enc$values[i]){
      endIdx <- endIdxs[i]
      startIdx <- endIdx - enc$lengths[i] + 1
      
      subx <- x[startIdx:endIdx]
      suby.max <- y.ci.max[startIdx:endIdx]
      suby.min <- y.ci.min[startIdx:endIdx]
      
      reduced.x <- c(subx, rev(subx))
      reduced.y <- c(suby.max, rev(suby.min))
      
      polygon(x=reduced.x, y=reduced.y, col=col, border=NA)
    }
  }
}

