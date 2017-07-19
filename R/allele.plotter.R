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
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. 
#' @param y.min.manual DEFAULT: NULL. Manually adds a min y-value.
#' @param no.title DEFAULT: FALSE. If TRUE, produces a plot with no title, not even the automatic
#' title components like formula and sample size.
#' @param override.title DEFAULT: NULL. If non-null, overrides automatic title components like
#' formula and sample size.
#' @param add.chr.to.label DEFAULT: FALSE. If TRUE, "chr" precedes integers on the x-axis.
#' @param alternative.labels DEFAULT: NULL. If non-null, specifies alternative labels for the
#' founder alleles.
#' @export
#' @examples allele.plotter.whole()
allele.plotter.whole <- function(scan.object, just.these.chr=NULL,
                                 scale="Mb", 
                                 main.colors=c(scales::alpha(rgb(240, 240, 0, maxColorValue=255), 0.8), # yellow
                                               scales::alpha(rgb(128, 128, 128, maxColorValue=255), 0.8), # gray
                                               scales::alpha(rgb(240, 128, 128, maxColorValue=255), 0.8), # pink/salmon
                                               scales::alpha(rgb(16, 16, 240, maxColorValue=255), 0.8), # dark blue
                                               scales::alpha(rgb(0, 160, 240, maxColorValue=255), 0.8), # light blue
                                               scales::alpha(rgb(0, 160, 0, maxColorValue=255), 0.8), # green
                                               scales::alpha(rgb(240, 0, 0, maxColorValue=255), 0.8), # red
                                               scales::alpha(rgb(144, 0, 224, maxColorValue=255), 0.8)), # purple,
                                 use.legend=TRUE, main="", my.bty="n", my.lwd=rep(1.5, 8),
                                 my.legend.cex=0.6, my.legend.pos="topright",
                                 y.max.manual=NULL, y.min.manual=NULL, no.title=FALSE, override.title=NULL,
                                 add.chr.to.label=FALSE, alternative.labels=NULL)
{
  allele.effects <- scan.object$allele.effects
  if(is.null(allele.effects)){
    stop("No allele effects in scan object. Re-run scan.h2lmm with return.allele.effects=TRUE", call.=FALSE)
  }
  chr <- scan.object$chr
  
  pos <- ifelse(rep(scale=="Mb", ncol(allele.effects)), scan.object$pos$Mb, scan.object$pos$cM)

  if(!is.null(just.these.chr)){
    keep.chr <- chr %in% just.these.chr
    chr <- chr[keep.chr]
    allele.effects <- allele.effects[, keep.chr]
    pos <- pos[keep.chr]
  }
  
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- max(as.numeric(unique(chr[chr != "X"]))) + 1
  }
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  allele.effects <- allele.effects[, order.i]
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x, na.rm=TRUE))
  max.pos <- tapply(pos, pre.chr, function(x) max(x, na.rm=TRUE))
  chr.types <- levels(pre.chr)
  
  # Finding max and min of y for plot window
  y.max <- ifelse(max(allele.effects) > 0, ceiling(max(allele.effects)), floor(max(allele.effects)))
  y.min <- ifelse(min(allele.effects) > 0, ceiling(min(allele.effects)), floor(min(allele.effects)))
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
  if(no.title){ this.title <- NULL }
  if(!is.null(override.title)){ this.title <- override.title }
  
  plot(pos[pre.chr==chr.types[1]], allele.effects[1, pre.chr==chr.types[1]], 
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)), 
       ylim=c(y.min, y.max), 
       xaxt="n", yaxt="n", xlab="", ylab="Additive allele effects", main=this.title,
       frame.plot=FALSE, type="l", pch=20, cex=0.5, lwd=1.5, col=main.colors[1])
  axis(side=2, at=y.min:y.max, las=2)
  
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos, na.rm=TRUE), min(this.pos, na.rm=TRUE):max(this.pos, na.rm=TRUE), max(this.pos, na.rm=TRUE)), 
                y=c(y.min, rep(y.max, length(min(this.pos, na.rm=TRUE):max(this.pos, na.rm=TRUE))), y.min), border=NA, col="gray88")
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)
      points(this.pos, allele.effects[1, pre.chr==chr.types[i]], type="l", lwd=1.5, col=main.colors[1])
      shift <- shift + max.pos[i]
    }
  }
  
  # Plot other alleles
  for(i in 2:nrow(allele.effects)){
    points(pos[pre.chr==chr.types[1]], allele.effects[i, pre.chr==chr.types[1]], type="l", col=main.colors[i], lwd=1.5)
    if(length(chr.types) > 1){
      shift <- max.pos[1]
      for(j in 2:length(chr.types)){
        points(pos[pre.chr==chr.types[j]] + shift, allele.effects[i, pre.chr==chr.types[j]], type="l", col=main.colors[i], lwd=1.5)
        shift <- shift + max.pos[j]
      }
    }
  }
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  if(!has.X){
    axis.label <- chr.types
  }
  if(add.chr.to.label){
    axis.label <- paste("chr", axis.label)
  }
  axis(side=1, tick=FALSE, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-1.5)
  if(use.legend){
    if(!is.null(alternative.labels)){ these.labels <- alternative.labels }
    else{ these.labels <- rownames(allele.effects) }
    legend(my.legend.pos, legend=these.labels, 
           lty=rep(1, nrow(allele.effects)), lwd=1.5, 
           col=main.colors, bty=my.bty, cex=my.legend.cex)
  }
}