#' Plot whole genome and single chromosome windows of a SNP-based genome scan
#'
#' This function takes the genome scan output from imputed.snp.scan.h2lmm() and plots the whole genome or single chromosome zoom-ins. 
#'
#' @param snp.scan An imputed.snp.scan.h2lmm() object.
#' @param just.these.chr DEFAULT: NULL. The chromosomes to be plotted. NULL leads to all chromosomes being plotted.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param title DEFAULT: "". Manually adds a max ylim to the plot. Allows multiple genome scans to easily be on the same scale.
#' @param alt.col DEFAULT: NULL. Allows for a custom color vector for individual SNPs.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param my.legend.pos DEFAULT: "topright". Specified position of the legend on the plot.
#' @param add.chr.to.label DEFAULT: FALSE. If TRUE, adds "Chr" before every chromosome label. If FALSE, "Chr" is added as an axis
#' label under the y-axis.
#' @param axis.cram DEFAULT: TRUE. This makes the plot much more likely to include all chromosome labels. With small plots, this could
#' lead to overlapping labels.
#' @param include.x.axis.line DEFAULT: TRUE. IF TRUE, this option adds an x-axis line with ticks between chromosomes.
#' @export
#' @examples snp.genome.plotter.whole()
snp.genome.plotter.whole <- function(snp.scan, just.these.chr=NULL, point.col="black", point.cex=0.5,
                                     scale="Mb", 
                                     distinguish.chr.type=c("color", "box"), distinguish.box.col="gray88", distinguish.snp.col="gray60",
                                     y.max.manual=NULL, my.y.line=2, my.y.axis.cex=1,
                                     title="", override.title=NULL, my.title.line=NA, title.cex=1,
                                     alt.col=NULL,
                                     hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL, thresholds.lty=2, thresholds.lwd=1,
                                     my.legend.cex=0.6, my.legend.pos="topright", my.bty="n", my.x.labels=TRUE,
                                     add.chr.to.label=FALSE, axis.cram=TRUE, include.x.axis.line=TRUE){
  
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  distinguish.chr.type <- distinguish.chr.type[1]
  main.object <- snp.scan
  outcome <- -log10(main.object$p.value)
  plot.this <- "p.value"
  this.ylab <- expression("-log"[10]*"P")
  
  chr <- main.object$chr
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
  # Allowing for special colors
  if(!is.null(alt.col)){ use.col <- alt.col }
  
  if(!is.null(just.these.chr)){
    keep.chr <- chr %in% just.these.chr
    chr <- chr[keep.chr]
    outcome <- outcome[keep.chr]
    pos <- pos[keep.chr]
  }
  
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  # Creating alternating colors for different chromosomes
  if(is.null(alt.col)){ 
    if(distinguish.chr.type == "color"){
      use.col <- c(point.col[1], distinguish.snp.col[1])[(as.numeric(droplevels(pre.chr)) %% 2 == 0) + 1]
    }
    else{
      use.col <- rep(point.col, length(outcome)) 
    }
  }
  
  outcome <- outcome[order.i]
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  use.col <- use.col[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x))
  max.pos <- tapply(pos, pre.chr, function(x) max(x))
  chr.types <- levels(pre.chr)
  
  # Finding max y of plot window
  y.max <- ceiling(max(outcome, hard.thresholds)) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  shift.left <- min(pos[chr==chr.types[1]])
  if(is.null(override.title)){
    condition.loci <- main.object$condition.loci
    if(!is.null(condition.loci)){
      this.title <- c(title,
                      paste(main.object$formula, paste(condition.loci, collapse=" + "),
                            paste0("SNP (", main.object$model.type, ")"), sep=" + "),
                      paste("n =", length(main.object$fit0$y)))
    }
    else{
      this.title <- c(title,
                      paste(main.object$formula, paste0("SNP (", main.object$model.type, ")"), sep=" + "),
                      paste("n =", length(main.object$fit0$y)))
    }
  }
  else{
    this.title <- override.title
  }
  
  x.max <- sum(max.pos)+(length(chr.types)-1)
  plot(1,
       xlim=c(0, x.max),
       ylim=c(-0.1, y.max),
       xaxt="n", yaxt="n", xlab="", ylab="", main=NA,
       frame.plot=FALSE, type="n")
  title(main=this.title, line=my.title.line, cex.main=title.cex)
  axis(side=2, at=0:y.max, las=1, cex.axis=my.y.axis.cex)
  mtext(text=this.ylab, side=2, line=my.y.line)
  
  label.spots <- max.pos[1]/2
  x.tick.spots <- c(0, max.pos[1])
  
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        if(distinguish.chr.type == "box"){
          polygon(x=c(min(this.pos), min(this.pos):max(this.pos), max(this.pos)), 
                  y=c(0, rep(y.max, length(min(this.pos):max(this.pos))), 0), border=NA, col=distinguish.box.col)
        }
      }
      label.spots <- c(label.spots, shift + max.pos[i]/2)
      x.tick.spots <- c(x.tick.spots, max.pos[i] + shift)
      
      shift <- shift + max.pos[i]
    }
  }
  # Adding in the points
  shift <- max.pos[1]
  
  if(length(chr.types) >= 1){
    points(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], pch=20, cex=point.cex, col=use.col[pre.chr==chr.types[1]])
    if(length(chr.types) > 1){
      for(i in 2:length(chr.types)){
        this.pos <- pos[pre.chr==chr.types[i]] + shift
        points(this.pos, outcome[pre.chr==chr.types[i]], type="p", pch=20, cex=point.cex, col=use.col[pre.chr==chr.types[i]])
        
        shift <- shift + max.pos[i]
      }
    }
  }
  
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  if(!has.X){
    axis.label <- chr.types
  }
  
  if(include.x.axis.line){
    axis(side=1, tick=TRUE, line=NA, at=x.tick.spots, 
         labels=NA, xpd=TRUE)
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
    
    if(!my.x.labels){ even.axis.label <- FALSE; odd.axis.label <- FALSE }
    
    axis(side=1, tick=FALSE, line=NA, at=odd.label.spots, labels=odd.axis.label, cex.axis=0.7, padj=-1.5, xpd=TRUE)
    axis(side=1, tick=FALSE, line=NA, at=even.label.spots, labels=even.axis.label, cex.axis=0.7, padj=-1.5, xpd=TRUE)
  }
  else{
    if(!my.x.labels){ axis.label <- FALSE }
    axis(side=1, tick=FALSE, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-1.5, xpd=TRUE)
  }
  
  if(!is.null(hard.thresholds)){
    if(length(thresholds.lty) == 1 & length(hard.thresholds) > 1){
      thresholds.lty <- rep(thresholds.lty, length(hard.thresholds))
    }
    if(length(thresholds.lwd) == 1 & length(hard.thresholds) > 1){
      thresholds.lwd <- rep(thresholds.lwd, length(hard.thresholds))
    }
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=thresholds.lty[i], 
             lwd=thresholds.lwd[i])
    }
  }
  if(!is.null(thresholds.legend)){
    legend(my.legend.pos, legend=thresholds.legend, col=thresholds.col, lwd=rep(thresholds.lwd, length(thresholds.legend)),
           lty=rep(thresholds.lty, length(thresholds.legend)),
           bty=my.bty, cex=my.legend.cex, bg="white")
  }
}

#' Plot a single chromosome window of a SNP-based genome scan overlayed with r^2 information
#'
#' This function takes the outputs from imputed.snp.scan.h2lmm() and pairwise.cor.snp.scan() and plots 
#' the SNP association for a single chromosome, overlayed with r^2 information for a specified SNP. 
#'
#' @param snp.scan An imputed.snp.scan.h2lmm() object.
#' @param r2.object A pairwise.cor.snp.scan() object.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param title DEFAULT: "". Manually adds a max ylim to the plot. Allows multiple genome scans to easily be on the same scale.
#' @param alt.col DEFAULT: NULL. Allows for a custom color vector for individual SNPs.
#' @param this.cex DEFAULT: 1. Allows for the adjustment of the cex value for the main plot.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param my.legend.pos DEFAULT: "topright". Specified position of the legend on the plot.
#' @param r2.bounds DEFAULT: NULL. If NULL, no interval is depicted on the plot. If set to a value in [0,1], will include interval
#' based on the given r2 on the plot.
#' @export
#' @examples snp.genome.plotter.w.r2()
snp.genome.plotter.w.r2 <- function(snp.scan, r2.object,
                                    scale="Mb", zoom.in=FALSE, zoom.in.by=0.1,
                                    y.max.manual=NULL, my.y.line=2, my.y.axis.cex=1, my.ylab.cex = 1,
                                    title="", override.title=NULL, my.title.line=NA, title.cex=1,
                                    alt.col=NULL, this.cex=1,
                                    hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                                    my.legend.cex=0.6, my.legend.pos="topleft", thresholds.lty=2, thresholds.lwd=1, my.bty="n", 
                                    r2.bounds=NULL, bounds.col="gray", high.color="red", low.color="blue", add.outline=FALSE,
                                    include.ramp=TRUE, ramp.cex=0.7){
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  main.object <- snp.scan
  
  outcome <- -log10(main.object$p.value)
  plot.this <- "p.value"
  this.ylab <- expression("-log"[10]*"P")
  
  # Allowing for special colors
  if(is.null(alt.col)){ use.col <- rep("black", length(outcome)) }
  if(!is.null(alt.col)){ use.col <- alt.col }
  
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
  point.locus <- r2.object$point.locus
  point.locus.outcome <- outcome[point.locus == main.object$loci]
  point.locus.pos <- pos[point.locus == main.object$loci]
  
  chr <- r2.object$chr
  outcome <- outcome[main.object$chr == chr]
  pos <- pos[main.object$chr == chr]
  
  order.i <- order(pos)
  
  outcome <- outcome[order.i]
  pos <- pos[order.i]
  use.col <- use.col[order.i]
  
  min.pos <- min(pos)
  max.pos <- max(pos)
  
  # Finding max y of plot window
  y.max <- ceiling(max(outcome, hard.thresholds)) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  this.title <- c(title,
                  paste0(main.object$formula, " + SNP (", main.object$model.type, ")"),
                  paste("n =", length(main.object$fit0$y)))
  if(!is.null(r2.bounds)){
    this.title <- c(this.title, 
                    paste("r2 interval level:", r2.bounds))
  }
  if(!is.null(override.title)){
    this.title <- override.title
  }
  this.xlab <- paste0("Chr ", chr, " Position (", scale, ")")
  
  high2low <- colorRampPalette(c(high.color, low.color))
  these.colors <- rev(high2low(1000))
  
  r2.col <- these.colors[ceiling(r2.object$r2*999.1)]
  
  if(!is.null(r2.bounds)){
    r2.interval <- extract.r2.interval(scan.object=snp.scan, r2.scan.object=r2.object, r2.level=r2.bounds)
    if(scale == "cM"){
      low.locus.pos <- r2.interval$lb.cM
      high.locus.pos <- r2.interval$ub.cM
    }
    else if(scale == "Mb"){
      low.locus.pos <- r2.interval$lb.Mb
      high.locus.pos <- r2.interval$ub.Mb
    }
    if(zoom.in){
      span <- high.locus.pos - low.locus.pos
      min.pos <- floor(low.locus.pos - zoom.in.by*span)
      max.pos <- ceiling(high.locus.pos + zoom.in.by*span)
      
      # Handling the boundaries
      min.pos <- max(min(pos), min.pos)
      max.pos <- min(max(pos), max.pos)
      
      new.span <- max.pos - min.pos
      
      cat(min.pos, "\n", max.pos, "\n")
    }
  }
  
  plot(0, pch='',
       xlim=c(min.pos, max.pos),
       ylim=c(0, y.max+1),
       yaxt="n", xlab=this.xlab, ylab="", main=NA,
       frame.plot=FALSE)
  title(main=this.title, line=my.title.line, cex.main=title.cex)
  axis(side=2, at=0:y.max, las=2, cex.axis=my.y.axis.cex)
  mtext(text=this.ylab, side=2, line=my.y.line, cex = my.ylab.cex)
  if(!is.null(r2.bounds)){
    polygon(c(rep(low.locus.pos, 2), rep(high.locus.pos, 2)), c(0, rep(y.max, 2), 0), col=bounds.col, border=NA)
  }
  if(add.outline){
    points(x=pos, y=outcome, col="black", bg=r2.col, pch=21, cex=this.cex)
  }
  else{
    points(x=pos, y=outcome, col=r2.col, pch=20, cex=this.cex)
  }
  
  points(x=point.locus.pos, y=point.locus.outcome, 
         bg=high.color, pch=21, cex=1.5)
  #axis(side=2, at=0:y.max, las=2)
  
  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=thresholds.lty, lwd=thresholds.lwd)
    }
  }
  if(!is.null(thresholds.legend)){
    legend(my.legend.pos, legend=thresholds.legend, col=thresholds.col, lty=rep(thresholds.lty, length(thresholds.legend)),
           lwd=rep(thresholds.lwd, length(thresholds.legend)), bty=my.bty, cex=my.legend.cex)
  }
  if(include.ramp){
    plotrix::color.legend(xl=max.pos - 0.25*new.span, yb=y.max, xr=max.pos, yt=y.max+0.5, cex=ramp.cex,
                          legend=c(0, 0.5, 1), rect.col=these.colors, align="rb", gradient="x")  
    text(x=max.pos - 0.125*new.span,
         y=y.max+0.75,
         labels="r2 with peak SNP")
  }
}