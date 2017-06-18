#' Plot whole genome and single chromosome windows of haplotype-based genome scan in a PDF output document
#'
#' This function takes the genome scan output from scan.h2lmm() and plots the whole genome and single chromosome zoom-ins for 
#' all the specified chromosomes. When multiple imputations are used, includes the 95\% confidence band on the median in the zoomed-in
#' plots.
#'
#' @param scan.object A scan.h2lmm() object (ROP or multiple imputations). If multiple imputations, median and confidence interval 
#' on median are plotted.
#' @param chr DEFAULT: c(1:19, "X"). The chromosomes to be plotted. DEFAULT is all the mouse chromosomes.
#' @param use.lod DEFAULT: FALSE. Plots either the LOD score or the -log10 p-value.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param median.band.col DEFAULT: "gray88". The color of the 95\% confident band plotted around the median.
#' @param main DEFAULT: "". Adds a title above the model.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' no threshols legend is used.
#' @param pdf.output.path That path of the PDF file to be generated.
#' @param pdf.height DEFAULT: 5. The height of an individual pages of the PDF.
#' @param pdf.width DEFAULT: 9. The width of an individual pages of the PDF.
#' @export
#' @examples genome.plotter.to.pdf()
genome.plotter.to.pdf <- function(scan.object, chr=c(1:19, "X"), use.lod=FALSE,
                                  scale=c("Mb", "cM"), main.col="black", median.band.col="gray88", main="", y.max.manual=NULL,
                                  hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                                  pdf.output.path, pdf.height=5, pdf.width=9, ...){
  scale <- scale[1]
  pdf(pdf.output.path, height=pdf.height, width=pdf.width)
  genome.plotter.whole(scan.list=list(scan.object), use.lod=use.lod,
                       scale=scale, main.colors=main.col, use.legend=FALSE,
                       main=main,
                       y.max.manual=y.max.manual,
                       hard.thresholds=hard.thresholds, thresholds.col=thresholds.col, thresholds.legend=thresholds.legend, ...)
  for(i in 1:length(chr)){
    genome.plotter.chr(scan.object=scan.object, chr=chr[i], use.lod=use.lod,
                       scale=scale, main.col=main.col, median.band.col=median.band.col,
                       main=main, y.max.manual=y.max.manual, 
                       hard.thresholds=hard.thresholds, thresholds.col=thresholds.col, thresholds.legend=thresholds.legend, ...)
  }
  dev.off()
}

#' Plot single chromosome windows of haplotype-based genome scan
#'
#' This function takes the genome scan output from scan.h2lmm() and plots the portion that corresponds to a single chromosome.
#' When multiple imputations are used, includes the 95\% confidence band on the median.
#'
#' @param scan.object A scan.h2lmm() object (ROP or multiple imputations). If multiple imputations, median and confidence interval 
#' on median are plotted.
#' @param chr The chromosome to be plotted.
#' @param use.lod DEFAULT: FALSE. Plots either the LOD score or the -log10 p-value.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM can be used.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param median.band.col DEFAULT: "gray88". The color of the 95\% confident band plotted around the median.
#' @param main DEFAULT: "". Adds a title above the model.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' no threshols legend is used.
#' @export
#' @examples genome.plotter.chr()
genome.plotter.chr <- function(scan.object, chr, use.lod=FALSE,
                               scale=c("Mb", "cM"), main.col="black", median.band.col="gray88",
                               main="",
                               y.max.manual=NULL, my.legend.cex=0.6,
                               hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL){
  scale <- scale[1]
  MI <- all.CI <- CI <- NULL
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  
  if(use.lod){
    all.outcome <- scan.object$LOD
    outcome <- all.outcome[scan.object$chr == chr]
    plot.this <- "LOD"
    this.ylab <- "LOD"
    if(!is.null(scan.object$MI.LOD)){
      all.MI <- scan.object$MI.LOD
      # Finding the 95% CI on the median
      all.CI <- apply(all.MI, 2, function(x) ci.median(x))
      CI <- all.CI[,scan.object$chr == chr]
    }
  }
  else{
    all.outcome <- -log10(scan.object$p.value)
    outcome <- all.outcome[scan.object$chr == chr]
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
    if(!is.null(scan.object$MI.LOD)){
      all.MI <- -log10(scan.object$MI.p.value)
      # Finding the 95% CI on the median
      all.CI <- apply(all.MI, 2, function(x) ci.median(x, conf=0.95))
      CI <- all.CI[,scan.object$chr == chr]
    }
  }
  
  if(scale == "Mb"){
    pos <- scan.object$pos$Mb[scan.object$chr == chr]
  }
  else if(scale == "c<"){
    pos <- scan.object$pos$cM[scan.object$chr == chr]
  }
  
  order.i <- order(pos)
  
  outcome <- outcome[order.i]
  pos <- pos[order.i]
  
  min.pos <- min(pos)
  max.pos <- max(pos)

  
  # Finding max y of plot window
  y.max <- ceiling(max(all.outcome, hard.thresholds, all.CI[2,])) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  this.title <- c(main, paste0(scan.object$formula, " + locus (", scan.object$model.type, ")"))
  
  plot(pos, outcome, 
       xlim=c(0, max.pos), 
       ylim=c(0, y.max), 
       yaxt="n", xlab=paste("Chr", chr, paste0("(", scale, ")")), ylab=this.ylab, main=this.title,
       frame.plot=F, type="l", lwd=1.5, col=main.col)
  axis(side=2, at=0:y.max, las=2)
  if(!is.null(CI)){
    polygon(x=c(pos, rev(pos)), y=c(CI[1,], rev(CI[2,])), density=NA, col=median.band.col)
  }
  points(pos, outcome, type="l", pch=20, cex=0.5, lwd=1.5, col=main.col)

  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2)
    }
  }
  if(!is.null(thresholds.legend)){
    legend("topleft", legend=thresholds.legend, col=thresholds.col, lty=rep(2, length(thresholds.legend)),
           bty="n", cex=my.legend.cex)
  }
}

#' Plot one or more haplotype-based genome scans flexibly (whole genomes or subset of chromosomes)
#'
#' This function takes the genome scan output from scan.h2lmm() and flexibly plots out the genome scan.
#'
#' @param scan.list A list of scan.h2lmm() objects that are to be plotted within the same genome scan plot.
#' @param use.lod DEFAULT: FALSE. Plots either the LOD score or the -log10 p-value.
#' @param just.these.chr DEFAULT: NULL. Specifies a subset of the chromosomes to be plotted. NULL results in all chromosomes being plotted.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param use.legend DEFAULT: TRUE. Include a legend for the different associations. If TRUE, the labels are the names of the non.mi.scan.list object.
#' @param main DEFAULT: NULL. Adds a title above the model.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param my.legend.lwd DEFAULT: NULL. If NULL, all lines have lwd=1.5. If not, option specifies the lwds.
#' @param my.legend.pos DEFAULT: "topright". Specifies where to put the legend, if specified in use.legend.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' no threshols legend is used.
#' @export
#' @examples genome.plotter.whole()
genome.plotter.whole <- function(scan.list, use.lod=FALSE, just.these.chr=NULL,
                                 scale="Mb", main.colors=c("black", "gray48", "blue"),
                                 use.legend=TRUE, main="",
                                 my.legend.cex=0.6, my.legend.lwd=NULL, my.legend.pos="topright",
                                 y.max.manual=NULL,
                                 hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL)
{
  if(is.null(my.legend.lwd)){ my.legend.lwd <- rep(1.5, length(scan.list)) }
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  main.object <- scan.list[[1]]
  if(use.lod){
    outcome <- main.object$LOD
    plot.this <- "LOD"
    this.ylab <- "LOD"
  }
  if(!use.lod){
    outcome <- -log10(main.object$p.value)
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
  }
  chr <- main.object$chr
  
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
  if(!is.null(just.these.chr)){
    keep.chr <- chr %in% just.these.chr
    chr <- chr[keep.chr]
    outcome <- outcome[keep.chr]
    pos <- pos[keep.chr]
  }
  
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- max(as.numeric(unique(chr[chr != "X"]))) + 1
  }
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  outcome <- outcome[order.i]
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x, na.rm=TRUE))
  max.pos <- tapply(pos, pre.chr, function(x) max(x, na.rm=TRUE))
  chr.types <- levels(pre.chr)
  
  # Finding max y of plot window
  y.max <- ceiling(max(outcome, hard.thresholds)) 
  if(length(scan.list) > 1){
    for(i in 2:length(scan.list)){
      if(use.lod){
        y.max <- ceiling(max(y.max, unlist(scan.list[[i]][plot.this])))
      }
      if(!use.lod){
        y.max <- ceiling(max(y.max, -log10(unlist(scan.list[[i]][plot.this]))))
      }
    }
  }
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  shift.left <- min(pos[chr==chr.types[1]], na.rm=TRUE)

  ### Handling the annoying differences between lmer and lm objects
  if(class(scan.list[[1]]$fit0) != "lmerMod"){
    this.title <- c(main, 
                    paste0(scan.list[[1]]$formula, " + locus (", scan.list[[1]]$model.type, ")"),
                    paste("n =", ifelse(is.null(scan.list[[1]]$fit0$weights), 
                                        length(scan.list[[1]]$fit0$y),
                                        sum(scan.list[[1]]$fit0$weights))))
  }
  else{
    this.title <- c(main, 
                    paste0(scan.list[[1]]$formula, " + locus (", scan.list[[1]]$model.type, ")"),
                    paste("n =", sum(scan.list[[1]]$fit0@resp$weights)))
  }
  
  plot(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], 
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)), 
       ylim=c(-0.1, y.max), 
       xaxt="n", yaxt="n", xlab="", ylab=this.ylab, main=this.title,
       frame.plot=F, type="l", pch=20, cex=0.5, lwd=my.legend.lwd[1], col=main.colors[1])
  axis(side=2, at=0:y.max, las=2)
  
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos, na.rm=TRUE), min(this.pos, na.rm=TRUE):max(this.pos, na.rm=TRUE), max(this.pos, na.rm=TRUE)), 
                y=c(0, rep(y.max, length(min(this.pos, na.rm=TRUE):max(this.pos, na.rm=TRUE))), 0), border=NA, col="gray88")
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)
      points(this.pos, outcome[pre.chr==chr.types[i]], type="l", lwd=my.legend.lwd[1], col=main.colors[1])
      shift <- shift + max.pos[i]
    }
  }
  
  # Plot other method's statistics
  if(length(scan.list) > 1){
      for(i in 2:length(scan.list)){
        
        this.scan <- scan.list[[i]]
        if(use.lod){
          compar.outcome <- this.scan$LOD
        }
        if(!use.lod){
          compare.outcome <- -log10(this.scan$p.value)
        }
        pos <- ifelse(rep(scale=="Mb", length(compare.outcome)), this.scan$pos$Mb, this.scan$pos$cM)
        
        ## Resetting for new scan objects
        chr <- this.scan$chr
        if(!is.null(just.these.chr)){
            keep.chr <- chr %in% just.these.chr
            chr <- chr[keep.chr]
            compare.outcome <- compare.outcome[keep.chr]
            pos <- pos[keep.chr]
        }
        
        has.X <- FALSE
        if(any(chr=="X")){
          has.X <- TRUE
          chr[chr=="X"] <- max(as.numeric(unique(chr[chr != "X"]))) + 1
        }
        
        pre.chr <- as.factor(as.numeric(chr))
        order.i <- order(pre.chr, pos)
        
        compare.outcome <- compare.outcome[order.i]
        pre.chr <- pre.chr[order.i]
        pos <- pos[order.i]
        
        min.pos <- tapply(pos, pre.chr, function(x) min(x, na.rm=TRUE))
        max.pos <- tapply(pos, pre.chr, function(x) max(x, na.rm=TRUE))
        chr.types <- levels(pre.chr)
        
        compare.shift <- max.pos[1]
        points(pos[pre.chr==chr.types[1]], compare.outcome[pre.chr==chr.types[1]], type="l", col=main.colors[i], lwd=my.legend.lwd[i])
        if(length(chr.types) > 1){
          for(j in 2:length(chr.types)){
            points(pos[pre.chr==chr.types[j]] + compare.shift, compare.outcome[pre.chr==chr.types[j]], type="l", col=main.colors[i], lwd=my.legend.lwd[i])
            compare.shift <- compare.shift + max.pos[j]
          }
        }
      }
    }
    if(has.X){
      axis.label <- c(chr.types[-length(chr.types)], "X")
    }
    if(!has.X){
      axis.label <- chr.types
    }
    axis(side=1, tick=F, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-1.5)
  if(use.legend){
    legend(my.legend.pos, legend=names(scan.list), 
           lty=rep(1, length(scan.list)), lwd=my.legend.lwd, 
           col=main.colors[1:length(scan.list)], bty="n", cex=my.legend.cex)
  }
  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2)
    }
  }
  if(!is.null(thresholds.legend)){
    legend("topleft", legend=thresholds.legend, col=thresholds.col, lty=rep(2, length(thresholds.legend)),
           bty="n", cex=my.legend.cex)
  }
}

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
#' @export
#' @examples snp.genome.plotter.whole()
snp.genome.plotter.whole <- function(snp.scan, just.these.chr=NULL,
                                     scale="Mb",
                                     y.max.manual=NULL, title="", alt.col=NULL,
                                     hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                                     my.legend.cex=0.6, my.legend.pos="topright"){
  
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  main.object <- snp.scan
  outcome <- -log10(main.object$p.value)
  plot.this <- "p.value"
  this.ylab <- expression("-log"[10]*"P")
  
  # Allowing for special colors
  if(is.null(alt.col)){ use.col <- rep("black", length(outcome)) }
  if(!is.null(alt.col)){ use.col <- alt.col }
  chr <- main.object$chr
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
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
  this.title <- c(title,
                  paste0(main.object$formula, " + SNP (", main.object$model.type, ")"),
                  paste("n =", length(main.object$fit0$y)))

  plot(1,
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)),
       ylim=c(-1, y.max),
       xaxt="n", yaxt="n", xlab="", ylab=this.ylab, main=this.title,
       frame.plot=F, type="n")
  axis(side=2, at=0:y.max, las=2)
  
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos), min(this.pos):max(this.pos), max(this.pos)), 
                y=c(0, rep(y.max, length(min(this.pos):max(this.pos))), 0), border=NA, col="gray88")
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)

      shift <- shift + max.pos[i]
    }
  }
  # Adding in the points
  shift <- max.pos[1]
  
  if(length(chr.types) >= 1){
    points(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], pch=20, cex=0.5, col=use.col[pre.chr==chr.types[1]])
    if(length(chr.types) > 1){
      for(i in 2:length(chr.types)){
        this.pos <- pos[pre.chr==chr.types[i]] + shift
        points(this.pos, outcome[pre.chr==chr.types[i]], type="p", pch=20, cex=0.5, col=use.col[pre.chr==chr.types[i]])
        
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
  axis(side=1, tick=F, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-3.5)
  
  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2)
    }
  }
  if(!is.null(thresholds.legend)){
    legend(my.legend.pos, legend=thresholds.legend, col=thresholds.col, lty=rep(2, length(thresholds.legend)),
           bty="n", cex=my.legend.cex)
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
                                    scale="Mb",
                                    y.max.manual=NULL, title="", alt.col=NULL, this.cex=1,
                                    hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                                    my.legend.cex=0.6, my.legend.pos="topleft",
                                    r2.bounds=NULL){
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
  this.xlab <- paste0("Chr ", chr, " Position (", scale, ")")
  
  red2blue <- colorRampPalette(c("red", "blue"))
  these.colors <- rev(red2blue(1000))
  
  r2.col <- these.colors[ceiling(r2.object$r2*999.1)]
  
  plot(0, pch='',
       xlim=c(min.pos, max.pos),
       ylim=c(0, y.max+1),
       yaxt="n", xlab=this.xlab, ylab=this.ylab, main=this.title,
       frame.plot=F)
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
    polygon(c(rep(low.locus.pos, 2), rep(high.locus.pos, 2)), c(0, rep(y.max, 2), 0), col="gray", border=NA)
  }
  points(x=pos, y=outcome, col=r2.col, pch=20, cex=this.cex)
  points(x=point.locus.pos, y=point.locus.outcome, 
         bg="red", pch=21, cex=1.5)
  axis(side=2, at=0:y.max, las=2)
  
  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2)
    }
  }
  if(!is.null(thresholds.legend)){
    legend(my.legend.pos, legend=thresholds.legend, col=thresholds.col, lty=rep(2, length(thresholds.legend)),
           bty="n", cex=my.legend.cex)
  }
  
  plotrix::color.legend(xl=floor(0.75*max.pos), yb=y.max, xr=max.pos, yt=y.max+0.5, legend=c(0, 0.25, 0.5, 0.75, 1), rect.col=these.colors, align="rb", gradient="x")  
  text(x=(max.pos - floor(0.75*max.pos))/2 + floor(0.75*max.pos),
       y=y.max+0.75,
       labels="r2 with peak SNP")
}

#' Plot whole genome and single chromosome windows of haplotype-based genome scan in a PDF output document
#'
#' This function takes the genome scan output from scan.h2lmm() and plots the whole genome and single chromosome zoom-ins for 
#' all the specified chromosomes. When multiple imputations are used, includes the 95\% confidence band on the median in the zoomed-in
#' plots.
#'
#' @param scan.object A scan.h2lmm() object (ROP or multiple imputations). If multiple imputations, median and confidence interval 
#' on median are plotted. Expected to the scan of the actual data.
#' @param qtl.ci.object A run.positional.scans() object. Should contain single chromosome scans from some form of sampling process, such as a parametric bootstrap.
#' @param ci.type Positional confidence interval to be included in the title. Example: "Parametric Bootstrap".
#' @param scan.type Scan type to be included in the title. Example: "ROP".
#' @param these.col Colors to be used for individual artificial scans.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either "Mb" or "cM" is expected.
#' @export
#' @examples single.chr.plotter.w.ci()
single.chr.plotter.w.ci <- function(scan.object, qtl.ci.object, 
                                    ci.type, scan.type, 
                                    these.col=c("#7BAFD4", "red"), scale="Mb"){
  
  outcome <- -log10(scan.object$p.value[scan.object$chr == qtl.ci.object$chr]) 
  this.ylab <- expression("-log"[10]*"P")
  
  pos <- qtl.ci.object$full.results$pos[[scale]]
  
  order.i <- order(pos)
  pos <- pos[order.i]
  outcome <- outcome[order.i]
  
  all.loci <- colnames(qtl.ci.object$full.results$p.values)
  
  peak.pos <- qtl.ci.object$peak.pos[[scale]]
  loci <- qtl.ci.object$peak.loci
  
  # Process per CI
  ci <- qtl.ci.object$ci[[scale]]

  lb.dist <- pos - ci[1]
  low.locus <- all.loci[lb.dist <= 0][which.max(lb.dist[lb.dist <= 0])]
  low.locus.pos <- pos[which(all.loci == low.locus)]
  ub.dist <- pos - ci[2]
  high.locus <- all.loci[ub.dist >= 0][which.min(ub.dist[ub.dist >= 0])]
  high.locus.pos <- pos[which(all.loci == high.locus)]
    
  actual.p.value <- scan.object$p.value[scan.object$chr == qtl.ci.object$chr]
  actual.loci <- scan.object$loci[scan.object$chr == qtl.ci.object$chr]
  actual.pos <- scan.object$pos[[scale]][scan.object$chr == qtl.ci.object$chr]
  region <- actual.pos >= low.locus.pos & actual.pos <= high.locus.pos
  peak.locus <- all.loci[region][which.min(actual.p.value[region])]
  peak.locus.pos <- actual.pos[region][which.min(actual.p.value[region])]
  
  main.title <- c(paste0(scan.type, ": ", scan.object$formula, " + locus (", scan.object$model.type, ")"),
                  paste0("QTL interval type: ", ci.type),
                  paste0("Width: ", round(high.locus.pos - low.locus.pos, 2), scale),
                  paste0("peak locus: ", peak.locus, " (", round(peak.locus.pos, 3), scale, ")"),
                  paste0("(closest) lower locus: ", low.locus, " (", round(low.locus.pos, 3), scale, ")"),
                  paste0("(closest) upper locus: ", high.locus, " (", round(high.locus.pos, 3), scale, ")"))
  full.results <- qtl.ci.object$full.results$p.values
  this.xlab <- paste("Chr", qtl.ci.object$chr, paste0("(", scale, ")"))
  if(!is.null(full.results)){ y.max <- max(outcome, -log10(full.results)) }
  else{ y.max <- max(outcome) }
  plot(1, 
       xlim=c(0, max(pos)), 
       ylim=c(0, y.max), 
       xlab=this.xlab, ylab=this.ylab, main=main.title,
       frame.plot=FALSE, type="l", pch=20, cex=0.5, las=1, cex.main=0.8)
  polygon(c(rep(low.locus.pos, 2), rep(high.locus.pos, 2)), c(0, rep(ceiling(max(outcome)), 2), 0), col="gray", border=NA)
  peaks <- qtl.ci.object$peak.loci.pos[[scale]]
  
  for(i in 1:nrow(full.results)){
    lines(pos, -log10(full.results[i,]), lwd=0.5, col=scales::alpha(these.col[i], 0.5))
  }
  rug(qtl.ci.object$peak.loci.pos[[scale]], col=scales::alpha("black", 0.5))

  lines(pos, outcome, lwd=1.5)
}

inspect.ci.genome.plotter.whole <- function(ci.object, scan.type.label, which.ci=1, ...){
  this.scan <- list(p.value = ci.object$full.results[which.ci,],
                    chr=rep(ci.object$chr, length(ci.object$full.results[which.ci,])),
                    pos = ci.object$pos)
  this.scan.list <- list()
  this.scan.list[[scan.type.label]] <- this.scan
  genome.plotter.whole(scan.list=this.scan.list, use.lod=FALSE, scale="cM", ...)
}


