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
#' @param no.title DEFAULT: FALSE. If TRUE, no title is printed.
#' @param override.title DEFAULT: NULL. If a string is specified, it is included on plot without any of the default automated title.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' no thresholds legend is used.
#' @param pdf.output.path That path of the PDF file to be generated.
#' @param pdf.height DEFAULT: 5. The height of an individual pages of the PDF.
#' @param pdf.width DEFAULT: 9. The width of an individual pages of the PDF.
#' @export
#' @examples genome.plotter.to.pdf()
genome.plotter.to.pdf <- function(scan.object, chr=c(1:19, "X"), use.lod=FALSE,
                                  scale=c("Mb", "cM"), main.col="black", median.band.col="gray88", 
                                  main="", no.title=FALSE, override.title=NULL,
                                  y.max.manual=NULL,
                                  hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                                  pdf.output.path, pdf.height=5, pdf.width=9, ...){
  scale <- scale[1]
  pdf(pdf.output.path, height=pdf.height, width=pdf.width)
  genome.plotter.whole(scan.list=list(scan.object), use.lod=use.lod,
                       scale=scale, main.colors=main.col, use.legend=FALSE,
                       main=main, no.title=no.title, override.title=override.title,
                       y.max.manual=y.max.manual,
                       hard.thresholds=hard.thresholds, thresholds.col=thresholds.col, thresholds.legend=thresholds.legend, ...)
  for(i in 1:length(chr)){
    genome.plotter.chr(scan.object=scan.object, chr=chr[i], use.lod=use.lod,
                       scale=scale, main.col=main.col, median.band.col=median.band.col,
                       main=main, no.title=no.title, override.title=override.title,
                       y.max.manual=y.max.manual, 
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
#' @param main.col DEFAULT: "black". The color of the main association score to be plotted.
#' @param median.band.col DEFAULT: "gray88". The color of the 95\% confident band plotted around the median.
#' @param main DEFAULT: "". Adds a title above the model.
#' @param no.title DEFAULT: FALSE. If TRUE, no title is printed.
#' @param override.title DEFAULT: NULL. If a string is specified, it is included on plot without any of the default automated title.
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
                               main="", no.title=FALSE, override.title=NULL, 
                               my.y.line=2, my.y.axis.cex=1, y.max.manual=NULL, make.y.axis.sparse=FALSE, my.ylab.cex=1,
                               my.x.line=2, my.x.axis.cex=1, my.xlab.cex=1, x.padj=1,
                               my.x.labels=TRUE, override.xlab=NULL, 
                               my.title.line=0.5, title.cex=1,
                               my.lwd=1.5,
                               override.ylab=NULL,
                               my.legend.cex=0.6, my.type="l", point.cex=0.5,
                               hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                               include.qtl.rug=FALSE, rug.pos=NULL, rug.col="gray50",
                               physical.dist.is.Mb=TRUE){
  scale <- scale[1]
  MI <- all.CI <- CI <- NULL
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  
  if(use.lod){
    all.outcome <- scan.object$LOD
    outcome <- all.outcome[scan.object$chr == chr & !is.na(scan.object$pos[[scale]])]
    plot.this <- "LOD"
    this.ylab <- "LOD"
    if(!is.null(scan.object$MI.LOD)){
      all.MI <- scan.object$MI.LOD
      # Finding the 95% CI on the median
      all.CI <- apply(all.MI, 2, function(x) ci.median(x))
      CI <- all.CI[,scan.object$chr == chr & !is.na(scan.object$pos[[scale]])]
    }
  }
  else{
    all.outcome <- -log10(scan.object$p.value)
    outcome <- all.outcome[scan.object$chr == chr & !is.na(scan.object$pos[[scale]])]
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
    if(!is.null(scan.object$MI.LOD)){
      all.MI <- -log10(scan.object$MI.p.value)
      # Finding the 95% CI on the median
      all.CI <- apply(all.MI, 2, function(x) ci.median(x, conf=0.95))
      CI <- all.CI[,scan.object$chr == chr & !is.na(scan.object$pos[[scale]])]
    }
  }
  
  pos <- scan.object$pos[[scale]][scan.object$chr == chr & !is.na(scan.object$pos[[scale]])]
  if(scale == "Mb" & !physical.dist.is.Mb){ # If for some reason the recorded Mb are actually bp
    pos <- pos/1000000
  }
  
  order.i <- order(pos)
  
  outcome <- outcome[order.i]
  pos <- pos[order.i]
  
  min.pos <- min(pos, na.rm=TRUE)
  max.pos <- max(pos, na.rm=TRUE)

  # Finding max y of plot window
  y.max <- ceiling(max(all.outcome, hard.thresholds, all.CI[2,])) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  if(!is.null(scan.object$locus.effect.type)){
    locus.effect.type <- ifelse(scan.object$locus.effect.type == "fixed", "fixef", "ranef")
    locus.term <- paste("locus", locus.effect.type, sep=".")
  }
  else{
    locus.term <- "locus"
  }
  this.title <- c(main, paste0(scan.object$formula, " + ", locus.term, " (", scan.object$model.type, ")"))
  if(no.title){ this.title <- NULL }
  if(!is.null(override.title)){ this.title <- override.title }
  
  plot(pos, outcome, 
       xlim=c(0, max.pos), 
       ylim=c(0, y.max), 
       xaxt="n", yaxt="n", xlab="", ylab="", main="",
       frame.plot=FALSE, type=my.type, cex=point.cex, lwd=my.lwd, col=main.col, pch=20)
  title(main=this.title, line=my.title.line, cex.main=title.cex)
  
  if(!make.y.axis.sparse){
    axis(side=2, at=0:y.max, las=2, cex.axis=my.y.axis.cex)
  }
  else{
    axis(side=2, las=2, cex.axis=my.y.axis.cex)
  }
  this.ylab <- ifelse(is.null(override.ylab), this.ylab, override.ylab)
  mtext(text=this.ylab, side=2, line=my.y.line, cex=my.ylab.cex)
  
  this.xlab <- ifelse(is.null(override.xlab), paste("Chr", chr, paste0("(", scale, ")")), override.xlab)
  axis(side=1, cex.axis=my.x.axis.cex, labels=my.x.labels, padj=x.padj)
  mtext(text=this.xlab, side=1, line=my.x.line, cex=my.xlab.cex)
  
  if(!is.null(CI)){
    polygon(x=c(pos, rev(pos)), y=c(CI[1,], rev(CI[2,])), density=NA, col=median.band.col)
  }
  points(pos, outcome, type=my.type, pch=20, cex=0.5, lwd=1.5, col=main.col)

  if(include.qtl.rug){
    if(is.null(rug.pos)){ rug.pos <- pos[which.max(outcome)] }
    if(length(rug.col) == 1){ rug.col <- rep(rug.col, length(rug.pos))}
    for(i in 1:length(rug.pos)){
      axis(1, at=rug.pos[i], col.ticks=rug.col[i], label=FALSE, cex.axis=0.1, lwd.ticks=3, lend="butt", tck=-0.1)
    }
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

#' Plot one or more haplotype-based genome scans flexibly (whole genomes or subset of chromosomes)
#'
#' This function takes the genome scan output from scan.h2lmm() and flexibly plots out the genome scan.
#'
#' @param scan.list A list of scan.h2lmm() objects that are to be plotted within the same genome scan plot.
#' @param use.lod DEFAULT: FALSE. Plots either the LOD score or the -log10 p-value.
#' @param just.these.chr DEFAULT: NULL. Specifies a subset of the chromosomes to be plotted. NULL results in all chromosomes being plotted.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param distinguish.chr.type DEFAULT: "color". The default specifies rectangular block backgrounds to distinguish adjacent chromosomes. The "color" option 
#' specifies alternating colors.
#' @param distinguish.box.col DEFAULT: "gray88". If distinguish.chr.type="box" is specified, this argument provides the color of the background rectangle.
#' @param distinguish.chr.col DEFAULT: "gray60". If distinguish.chr.type="color" is specified, this argument provides the alternating color. Multiple
#' colors should be specified if multiple scans are included in the plot.
#' @param use.legend DEFAULT: TRUE. Include a legend for the different associations. If TRUE, the labels are the names of the non.mi.scan.list object.
#' @param main DEFAULT: NULL. Adds a title above the model.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param my.legend.lwd DEFAULT: NULL. If NULL, all lines have lwd=1.5. If not, option specifies the lwds.
#' @param my.legend.pos DEFAULT: "topright". Specifies where to put the legend, if specified in use.legend.
#' @param no.title DEFAULT: FALSE. If TRUE, no title is printed.
#' @param override.title DEFAULT: NULL. If a string is specified, it is included on plot without any of the default automated title.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' no threshols legend is used.
#' @param add.chr.to.label DEFAULT: FALSE. If TRUE, adds "Chr" before every chromosome label. If FALSE, "Chr" is added as an axis
#' label under the y-axis.
#' @param axis.cram DEFAULT: TRUE. This makes the plot much more likely to include all chromosome labels. With small plots, this could
#' lead to overlapping labels.
#' @param include.x.axis.line DEFAULT: TRUE. IF TRUE, this option adds an x-axis line with ticks between chromosomes.
#' @export
#' @examples genome.plotter.whole()
genome.plotter.whole <- function(scan.list, use.lod=FALSE, just.these.chr=NULL,
                                 scale="Mb", main.colors=c("black", "cyan", "darkgreen"),
                                 distinguish.chr.type=c("color", "box"), distinguish.box.col="gray88", 
                                 distinguish.chr.col=c("gray60", "#008080", "greenyellow"),
                                 use.legend=TRUE, main="",
                                 my.legend.cex=0.6, my.legend.lwd=NULL, my.legend.lty=1, my.legend.pos="topright",
                                 y.max.manual=NULL, my.y.line=2, my.y.axis.cex=1, my.y.lab.cex=0.7,
                                 my.x.lab.cex=0.7, my.x.labels=TRUE,
                                 no.title=FALSE, override.title=NULL, my.title.line=NA, title.cex=1,
                                 hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL,
                                 thresholds.lwd=NULL, thresholds.legend.pos="topleft", thresholds.lty=NULL,
                                 add.chr.to.label=FALSE, axis.cram=TRUE, include.x.axis.line=TRUE,
                                 mark.locus=NULL, mark.locus.col="red", which.mark=1,
                                 add.polygon=FALSE, which.polygon=1){
  # If list has no names, use.legend is set to FALSE
  if (is.null(names(scan.list))){ use.legend=FALSE }
  if (is.null(my.legend.lwd)){ my.legend.lwd <- rep(1.5, length(scan.list)) }
  if (is.null(thresholds.lwd)){ thresholds.lwd <- rep(1, ifelse(is.null(hard.thresholds), 0, length(hard.thresholds))) }
  if (is.null(thresholds.lty)){ thresholds.lty <- rep(2, ifelse(is.null(hard.thresholds), 0, length(hard.thresholds))) }
  if (length(my.legend.lty) == 1) { my.legend.lty <- rep(my.legend.lty, length(scan.list)) }
  
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  
  distinguish.chr.type <- distinguish.chr.type[1]
  main.object <- scan.list[[1]]
  if(use.lod){
    outcome <- main.object$LOD
    plot.this <- this.ylab <- "LOD"
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
  
  ### Fixef or ranef
  if(length(scan.list) == 1 & !is.null(scan.list[[1]]$locus.effect.type)){
    locus.effect.type <- ifelse(scan.list[[1]]$locus.effect.type == "fixed", "fixef", "ranef")
    locus.term <- paste("locus", locus.effect.type, sep=".")
  }
  else{
    locus.term <- "locus"
  }

  if(no.title){ this.title <- NULL }
  else if(!is.null(override.title)){ this.title <- override.title }
  else{
    ### Handling the annoying differences between lmer and lm objects
    if(is.null(scan.list[[1]]$fit0)){
      this.title <- c(main, 
                      paste0(scan.list[[1]]$formula, " + ", locus.term, " (", scan.list[[1]]$model.type, ")"),
                      paste("n =", round(scan.list[[1]]$n, 2)))
    }
    else{
      if(class(scan.list[[1]]$fit0) != "lmerMod"){
        this.title <- c(main, 
                        paste0(scan.list[[1]]$formula, " + ", locus.term, " (", scan.list[[1]]$model.type, ")"),
                        paste("n =", round(ifelse(is.null(scan.list[[1]]$fit0$weights), 
                                                  length(scan.list[[1]]$fit0$y),
                                                  sum(scan.list[[1]]$fit0$weights)), 2)))
      }
      else{
        this.title <- c(main, 
                        paste0(scan.list[[1]]$formula, " + ", locus.term, " (", scan.list[[1]]$model.type, ")"),
                        paste("n =", round(sum(scan.list[[1]]$fit0@resp$weights), 2)))
      }
    }
  }
  
  ## For later plotting, like mark.locus
  updated.pos <- rep(NA, length(outcome))
  names(updated.pos) <- names(outcome)
  updated.pos[pre.chr==chr.types[1]] <- pos[pre.chr==chr.types[1]]
  
  x.max <- sum(max.pos)+(length(chr.types)-1)
  plot(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], 
       xlim=c(0, x.max), 
       ylim=c(-0.1, y.max), 
       xaxt="n", yaxt="n", ylab="", xlab="", main=NA,
       frame.plot=FALSE, type="l", pch=20, cex=0.5, lwd=my.legend.lwd[1], lty=my.legend.lty[1], col=main.colors[1])
  title(main=this.title, line=my.title.line, cex.main=title.cex)
  axis(side=2, at=0:y.max, las=2, cex.axis=my.y.axis.cex)
  mtext(text=this.ylab, side=2, line=my.y.line, cex=my.y.lab.cex)
  
  if (1 %in% which.polygon & add.polygon) {
    polygon(expand.for.polygon.x(pos[pre.chr==chr.types[1]]), 
            expand.for.polygon.y(outcome[pre.chr==chr.types[1]]), 
            col=main.colors[1],
            border=NA)
  }
  
  label.spots <- max.pos[1]/2
  x.tick.spots <- c(0, max.pos[1])
  shift <- max.pos[1]
  
  if(distinguish.chr.type == "box"){ 
    this.col <- rep(main.colors[1], length(outcome)) 
  }
  else if(distinguish.chr.type == "color"){ 
    this.col <- c(main.colors[1], distinguish.chr.col[1])[(as.numeric(as.character(pre.chr)) %% 2 == 0) + 1]
  }
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      updated.pos[pre.chr==chr.types[i]] <- pos[pre.chr==chr.types[i]] + shift
      if(distinguish.chr.type == "box"){
        if(i %% 2 == 0){
          polygon(x=c(shift, max(this.pos, na.rm=TRUE), max(this.pos, na.rm=TRUE), shift), 
                  y=c(y.max, y.max, 0, 0), border=NA, col=distinguish.box.col)
          
        }
      }
      label.spots <- c(label.spots, shift + max.pos[i]/2)
      x.tick.spots <- c(x.tick.spots, max.pos[i] + shift)
      points(this.pos, outcome[pre.chr==chr.types[i]], type="l", lty=my.legend.lty[1], lwd=my.legend.lwd[1], col=this.col[pre.chr==chr.types[i]])
      if (1 %in% which.polygon & add.polygon) {
        polygon(expand.for.polygon.x(this.pos), 
                expand.for.polygon.y(outcome[pre.chr==chr.types[i]]), 
                col=this.col[pre.chr==chr.types[i]],
                border=NA)
      }
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
      
      ## Setting up colors for additional scans
      if(distinguish.chr.type == "box"){ 
        this.col <- rep(main.colors[i], length(compare.outcome)) 
      }
      else if(distinguish.chr.type == "color"){ 
        this.col <- c(main.colors[i], distinguish.chr.col[i])[(as.numeric(as.character(pre.chr)) %% 2 == 0) + 1]
      }
      
      compare.shift <- max.pos[1]
      points(pos[pre.chr==chr.types[1]], compare.outcome[pre.chr==chr.types[1]], type="l", 
             col=this.col[pre.chr==chr.types[1]], lwd=my.legend.lwd[i], lty=my.legend.lty[i])
      
      if (i %in% which.polygon & add.polygon) {
        polygon(expand.for.polygon.x(pos[pre.chr==chr.types[1]]), 
                expand.for.polygon.y(compare.outcome[pre.chr==chr.types[1]]), 
                col=this.col[pre.chr==chr.types[1]],
                border=NA)
      }
      
      ## For later plotting, like mark.locus
      if (i == which.mark) {
        updated.pos <- rep(NA, length(compare.outcome))
        names(updated.pos) <- names(compare.outcome)
        updated.pos[pre.chr==chr.types[1]] <- pos[pre.chr==chr.types[1]]
      }
      
      if (length(chr.types) > 1) {
        for (j in 2:length(chr.types)) {
          points(pos[pre.chr==chr.types[j]] + compare.shift, compare.outcome[pre.chr==chr.types[j]], type="l", 
                 col=this.col[pre.chr==chr.types[j]], 
                 lwd=my.legend.lwd[i],
                 lty=my.legend.lty[i])
          if (i %in% which.polygon & add.polygon) {
            polygon(expand.for.polygon.x(pos[pre.chr==chr.types[j]] + compare.shift), 
                    expand.for.polygon.y(compare.outcome[pre.chr==chr.types[j]]), 
                    col=this.col[pre.chr==chr.types[j]],
                    border=NA)
          }
          if (i == which.mark) {
            updated.pos[pre.chr==chr.types[j]] <- pos[pre.chr==chr.types[j]] + compare.shift
          }
          compare.shift <- compare.shift + max.pos[j]
        }
      }
    }
  }
  if (has.X) {
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  else {
    axis.label <- chr.types
  }
  
  if (include.x.axis.line) {
    axis(side=1, tick=TRUE, line=NA, at=x.tick.spots, 
         labels=NA, xpd=TRUE)
  }
  
  if (add.chr.to.label) {
    axis.label <- paste("Chr", axis.label)
  }
  else {
    axis.label <- c("Chr", axis.label)
    label.spots <- c(-0.04*x.max, label.spots)
  }
  
  if (axis.cram) {
    odd.axis.label <- axis.label[(1:length(axis.label) %% 2) == 1]
    odd.label.spots <- label.spots[(1:length(label.spots) %% 2) == 1]
    
    even.axis.label <- axis.label[(1:length(axis.label) %% 2) == 0]
    even.label.spots <- label.spots[(1:length(label.spots) %% 2) == 0]
    
    if (!my.x.labels) { even.axis.label <- FALSE; odd.axis.label <- FALSE }
    
    axis(side=1, tick=FALSE, line=NA, at=odd.label.spots, labels=odd.axis.label, cex.axis=my.x.lab.cex, padj=-1.5, xpd=TRUE)
    axis(side=1, tick=FALSE, line=NA, at=even.label.spots, labels=even.axis.label, cex.axis=my.x.lab.cex, padj=-1.5, xpd=TRUE)
  }
  else {
    if (!my.x.labels) { axis.label <- FALSE }
    axis(side=1, tick=FALSE, line=NA, at=label.spots, labels=axis.label, cex.axis=my.x.lab.cex, padj=-1.5, xpd=TRUE)
  }
  if (!is.null(mark.locus)) {
    rug(x=updated.pos[which(names(updated.pos) == mark.locus)], lwd=4, col=mark.locus.col)
  }
  if (use.legend) {
    if (add.polygon) {
      these.lty <- my.legend.lty
      these.lwd <- my.legend.lwd
      these.fill.col <- main.colors
      these.fill.border <- main.colors
      these.lty[which.polygon] <- these.lwd[which.polygon] <- NA
      these.fill.col[!(1:length(scan.list) %in% which.polygon)] <- these.fill.border[!(1:length(scan.list) %in% which.polygon)] <-  NA
      this.x.intersp=ifelse(is.na(these.fill.col), 2, 0.5)
      
      legend(my.legend.pos, 
             legend=names(scan.list), 
             lty=these.lty, 
             lwd=these.lwd,
             fill=these.fill.col,
             border=these.fill.border,
             x.intersp=this.x.intersp,
             col=main.colors[1:length(scan.list)], bty="n", cex=my.legend.cex)
    }
    else {
      legend(my.legend.pos, 
             legend=names(scan.list), 
             lty=my.legend.lty, 
             lwd=my.legend.lwd, 
             col=main.colors[1:length(scan.list)], bty="n", cex=my.legend.cex)
    }
  }
  if (!is.null(hard.thresholds)) {
    for (i in 1:length(hard.thresholds)) {
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=thresholds.lty[i], lwd=thresholds.lwd[i])
    }
  }
  if (!is.null(thresholds.legend)) {
    legend(thresholds.legend.pos, legend=thresholds.legend, col=thresholds.col, lty=thresholds.lty,
           lwd=thresholds.lwd, bty="n", cex=my.legend.cex)
  }
}

expand.for.polygon.x <- function(x){
  return(c(x[1], x, x[length(x)]))
}
expand.for.polygon.y <- function(x){
  return(c(0, x, 0))
}


#' Plot user-specified windows of haplotype-based and snp-based genome scans
#'
#' This function takes genome scan association outputs and plots the portion that corresponds to a region of a single chromosome.
#' When multiple imputations are used, includes the 95\% confidence band on the median.
#'
#' @param haplotype.association DEFAULT: NULL. A list of scan.h2lmm() objects (ROP or multiple imputations). If multiple imputations, median and confidence interval 
#' on median are plotted. If NULL, presumably only SNP scans will be plotted.
#' @param snp.association DEFAULT: NULL. A list of imputed.snp.scan.h2lmm() objects. If NULL, presumably only haplotype-based scans will be plotted.
#' @param use.lod DEFAULT: FALSE. Plots either the LOD score or the -log10 p-value.
#' @param chr The chromosome to be plotted.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM can be used.
#' @param region.min DEFAULT: NULL. The lower bound of the region to be plotted. Should match the scale. If NULL, defaults to the minimum of the specified chromosome.
#' @param region.max DEFAULT: NULL. The upper bound of the region to be plotted. Should match the scale. If NULL, defaults to the maximum of the specified chromosome.
#' @param haplotype.col DEFAULT: "black". The color of the haplotype-based association score to be plotted.
#' @param snp.col DEFAULT: "black". The color of the SNP association score to be plotted.
#' @param median.band.col DEFAULT: "gray88". The color of the 95\% confident band plotted around the median.
#' @param main DEFAULT: "". Adds a title above the model.
#' @param no.title DEFAULT: FALSE. If TRUE, no title is printed.
#' @param override.title DEFAULT: NULL. If a string is specified, it is included on plot without any of the default automated title.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param thresholds.legend DEFAULT: NULL. If non-NULL, string arguments used as labels in thresholds legend. If NULL,
#' no threshols legend is used.
#' @param use.legend DEFAULT: TRUE. Include a legend for the different associations. If TRUE, the labels are the names of the non.mi.scan.list object.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param my.legend.pos DEFAULT: "topright". Specifies where to put the legend, if specified in use.legend.
#' @export
#' @examples genome.plotter.region()
genome.plotter.region <- function(haplotype.association=NULL, snp.association=NULL, use.lod=FALSE,
                                  chr, region.min=NULL, region.max=NULL, scale=c("Mb", "cM"),
                                  haplotype.col=c("blue", "red"), 
                                  haplotype.lwd=3, haplotype.lty=1,
                                  median.band.col=c("cyan", "pink"),
                                  snp.col=c("black", "gray"), snp.pch=20, snp.cex=0.9,
                                  main="", no.title=FALSE, override.title=NULL,
                                  y.max.manual=NULL,
                                  my.y.line=2, my.y.axis.cex=1, my.y.lab.cex=0.5,
                                  my.x.line=2, my.x.axis.cex=1, my.xlab.cex=1, x.padj=-0.3,
                                  my.x.labels=TRUE, override.xlab=NULL, 
                                  my.title.line=0.5, my.title.cex=1,
                                  hard.thresholds=NULL, thresholds.col="red", thresholds.legend=NULL, thresholds.lwd=2,
                                  use.legend=TRUE, my.legend.cex=0.6, my.legend.pos="topright", my.bty="n",
                                  rug.pos=NULL, rug.col="gray50"){
  scale <- scale[1]

  if(is.null(haplotype.association) & is.null(snp.association)){
    stop("No associations included in arguments to genome.plotter.region()", call.=FALSE)
  }
  if(length(thresholds.col) < length(hard.thresholds)){ thresholds.col <- rep(thresholds.col, length(hard.thresholds)) }
  
  # Repeating parameters for multiple scans
  if(length(snp.pch) == 1 & length(snp.association) > 1){ 
    snp.pch <- rep(snp.pch, length(snp.association))
  }
  if(length(snp.cex) == 1 & length(snp.association) > 1){ 
    snp.cex <- rep(snp.cex, length(snp.association))
  }
  if(length(haplotype.lwd) == 1 & length(haplotype.association) > 1){ 
    haplotype.lwd <- rep(haplotype.lwd, length(haplotype.association))
  }
  if(length(haplotype.lty) == 1 & length(haplotype.association) > 1){
    haplotype.lty <- rep(haplotype.lty, length(haplotype.association))
  }
  if(length(thresholds.lwd) == 1 & length(hard.thresholds) > 1){
    thresholds.lwd <- rep(thresholds.lwd, length(hard.thresholds))
  }
  
  this.ylab <- ifelse(use.lod, "LOD", expression("-log"[10]*"P"))
  outcome.type <- ifelse(use.lod, "LOD", "p.value")
  
  ## Finding plot bounds and processing scans
  y.max <- x.min <- x.max <- haps.to.plot <- snps.to.plot <- CI <- NULL
  if(!is.null(haplotype.association)){
    haps.to.plot <- list()
    for(i in 1:length(haplotype.association)){
      this.scan <- haplotype.association[[i]]
      this.pos <- this.scan$pos[[scale]]
      this.outcome <- this.scan[[outcome.type]]
      if(!use.lod){ this.outcome <- -log10(this.outcome) }
      
      ## Getting index to plot
      keep.chr <- this.scan$chr == chr
      keep.na <- !is.na(this.pos)
      keep <- (keep.chr + keep.na) == 2
      order.i <- order(this.pos[keep])
      x.min <- min(x.min, 
                   grab.min.pos.from.scan(scan.object=this.scan, scale=scale, chr=chr))
      x.max <- max(x.max, 
                   grab.max.pos.from.scan(scan.object=this.scan, scale=scale, chr=chr))
      ## Imputations interval
      if(!is.null(this.scan$MI.LOD)){
        if(use.lod){
          all.MI <- this.scan$MI.LOD[,keep][,order.i]
        }
        else{
          all.MI <- -log10(this.scan$MI.p.value[,keep][,order.i])
        }
        # Finding the 95% CI on the median
        CI <- apply(all.MI, 2, function(x) ci.median(x, conf=0.95))
      }
      haps.to.plot[[i]] <- list(pos=this.pos[keep][order.i], 
                                outcome=this.outcome[keep][order.i],
                                CI=CI) 
    }
  }
  if(!is.null(snp.association)){
    snps.to.plot <- list()
    for(i in 1:length(snp.association)){
      this.scan <- snp.association[[i]]
      this.pos <- this.scan$pos[[scale]]
      this.outcome <- this.scan[[outcome.type]]
      keep.chr <- this.scan$chr == chr
      keep.na <- !is.na(this.pos)
      keep <- (keep.chr + keep.na) == 2
      order.i <- order(this.pos[keep])
      if(!use.lod){ this.outcome <- -log10(this.outcome) }
      x.min <- min(x.min, 
                   grab.min.pos.from.scan(scan.object=this.scan, scale=scale, chr=chr))
      x.max <- max(x.max, 
                   grab.max.pos.from.scan(scan.object=this.scan, scale=scale, chr=chr))
      snps.to.plot[[i]] <- list(pos=this.pos[keep][order.i], 
                                outcome=this.outcome[keep][order.i]) 
    }
  }
  if(!is.null(region.min)){ x.min <- region.min }
  if(!is.null(region.max)){ x.max <- region.max }
  
  ## Setting y limit within region
  if(!is.null(haplotype.association)){
    for(i in 1:length(haps.to.plot)){
      this.scan <- haps.to.plot[[i]]
      y.max <- max(y.max, max(this.scan$outcome[this.scan$pos >= x.min & this.scan$pos <= x.max], na.rm=TRUE), na.rm=TRUE)
      if(!is.null(haps.to.plot[[i]]$CI)){
        y.max <- max(y.max, max(haps.to.plot[[i]]$CI[,this.scan$pos >= x.min & this.scan$pos <= x.max], na.rm=TRUE), na.rm=TRUE)
      }
    }
  }
  if(!is.null(snp.association)){
    for(i in 1:length(snps.to.plot)){
      this.scan <- snps.to.plot[[i]]
      y.max <- max(y.max, max(this.scan$outcome[this.scan$pos >= x.min & this.scan$pos <= x.max], na.rm=TRUE), na.rm=TRUE)
    }
  }
  y.max <- max(y.max, max(c(hard.thresholds, 0)), na.rm=TRUE)
  if(!is.null(y.max.manual)){ y.max <- y.max.manual }
  
  ## Handling the title
  if(!is.null(haplotype.association[[1]])){ first.scan <- haplotype.association[[1]] }
  else{ first.scan <- snp.association[[1]] }
  if(!is.null(first.scan$locus.effect.type)){
    locus.effect.type <- ifelse(first.scan$locus.effect.type == "fixed", "fixef", "ranef")
    locus.term <- paste("locus", locus.effect.type, sep=".")
  }
  else{
    locus.term <- "locus"
  }
  this.title <- c(main, paste0(first.scan$formula, " + ", locus.term, " (", first.scan$model.type, ")"))
  if(no.title){ this.title <- NULL }
  if(!is.null(override.title)){ this.title <- override.title }
  
  plot(0, 
       xlim=c(x.min, x.max), 
       ylim=c(0, y.max), xaxt="n", yaxt="n",
       xlab="", ylab="", main="",
       frame.plot=FALSE, type="l", pch="", cex=0.5)
  title(main=this.title, line=my.title.line, cex.main=my.title.cex)
  
  x.ticks <- seq(x.min, x.max, length.out=5)
  x.ticks <- round(x.ticks)
  this.xlab <- ifelse(is.null(override.xlab), paste("Chr", chr, paste0("(", scale, ")")), override.xlab)
  axis(side=1, tick=TRUE, cex.axis=my.x.axis.cex, labels=my.x.labels, padj=x.padj)
  mtext(text=this.xlab, side=1, line=my.x.line, cex=my.xlab.cex)

  axis(side=2, las=2, cex.axis=my.y.axis.cex)
  mtext(text=this.ylab, side=2, line=my.y.line, cex=my.y.lab.cex)
  
  ## Adding associations
  if(!is.null(haplotype.association)){
    ## Plotting haplotype intervals
    for(i in 1:length(haplotype.association)){
      this.scan <- haps.to.plot[[i]]
      this.pos <- this.scan$pos
      
      this.col <- 1
      if(!is.null(this.scan$CI)){
        CI <- this.scan$CI
        polygon(x=c(this.pos, rev(this.pos)), y=c(CI[1,], rev(CI[2,])), density=NA, col=median.band.col[this.col])
      }
      else{
        this.col <- this.col + 1
      }
    }
    ## Plotting haplotype association lines
    for(i in 1:length(haplotype.association)){
      this.scan <- haps.to.plot[[i]]
      this.pos <- this.scan$pos
      this.outcome <- this.scan$outcome
      
      points(this.pos, this.outcome, col=haplotype.col[i], type="l", lty=haplotype.lty[i], lwd=haplotype.lwd[i])
    }
  }
  if(!is.null(snp.association)){
    for(i in 1:length(snp.association)){
      this.scan <- snps.to.plot[[i]]
      this.pos <- this.scan$pos
      this.outcome <- this.scan$outcome
      
      points(this.pos, this.outcome, 
             col=snp.col[i], pch=snp.pch[i], cex=snp.cex[i])
    }
  }
  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2, lwd=thresholds.lwd[i])
    }
  }
  if(!is.null(thresholds.legend)){
    legend("topleft", legend=thresholds.legend, col=thresholds.col, lty=rep(2, length(thresholds.legend)),
           bty="n", cex=my.legend.cex)
  }
  
  scan.names <- c(names(snp.association), names(haplotype.association))
  use.legend <- ifelse(is.null(scan.names), FALSE, TRUE)
  if(use.legend){
    if(is.null(snp.association)){ use.snp.col <- NULL }
    else{ use.snp.col <- snp.col[1:length(snp.association)] }
    if(is.null(haplotype.association)){ use.haplotype.col <- NULL }
    else{ use.haplotype.col <- haplotype.col[1:length(haplotype.association)] }
    legend(my.legend.pos, legend=scan.names, 
           lty=c(rep(NA, length(snp.association)), rep(haplotype.lty, length(haplotype.association))), 
           lwd=c(rep(NA, length(snp.association)), haplotype.lwd), 
           pch=c(rep(20, length(snp.association)), rep(NA, length(haplotype.association))),
           col=c(use.snp.col, use.haplotype.col), 
           bty=my.bty, bg="white", cex=my.legend.cex)
  }
  if(!is.null(rug.pos)){
    if(length(rug.col) == 1){ rug.col <- rep(rug.col, length(rug.pos))}
    for(i in 1:length(rug.pos)){
      axis(1, at=rug.pos[i], col.ticks=rug.col[i], label=FALSE, cex.axis=0.1, lwd.ticks=3, lend="butt", tck=-0.1)
    }
  }
}

grab.min.pos.from.scan <- function(scan.object, scale="Mb", chr=NULL){
  pos <- scan.object$pos[[scale]]
  chr.vec <- scan.object$chr
  if(!is.null(chr)){ pos <- pos[chr.vec %in% chr]}
  min.pos <- min(pos, na.rm=TRUE)
  return(min.pos)
}
grab.max.pos.from.scan <- function(scan.object, scale="Mb", chr=NULL){
  pos <- scan.object$pos[[scale]]
  chr.vec <- scan.object$chr
  if(!is.null(chr)){ pos <- pos[chr.vec %in% chr]}
  max.pos <- max(pos, na.rm=TRUE)
  return(max.pos)
}
grab.max.statistic.from.scan <- function(scan.object, outcome="p.value", chr=NULL){
  assoc <- scan.object[[outcome]]
  if(outcome == "p.value"){ assoc <- -log10(assoc)}
  chr.vec <- scan.object$chr
  if(!is.null(chr)){ assoc <- assoc[chr.vec %in% chr]}
  max.assoc <- max(assoc, na.rm=TRUE)
  return(max.assoc)
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
#' @param alpha DEFAULT: 0.05. The specified alpha level of the positional confidence interval.
#' @export
#' @examples single.chr.plotter.w.ci()
single.chr.plotter.w.ci <- function(scan.object, qtl.ci.object, 
                                    ci.type, scan.type, 
                                    these.col=c("#7BAFD4", "red"), scale="Mb",
                                    alpha=0.05){
  
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
  ci <- quantile(qtl.ci.object$peak.loci.pos[[scale]], probs=c(alpha/2, 1 - alpha/2))
  
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
                  paste0("QTL interval type: ", (1-alpha)*100, "% ", ci.type),
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

