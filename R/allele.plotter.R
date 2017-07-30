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
                                 use.legend=TRUE, main="", my.bty="n", my.lwd=rep(1.25, 8),
                                 set.plot.limit=c(-10, 10), # Null places no limit on y-axis
                                 my.legend.cex=0.6, my.legend.pos="topright", transparency=0.6,
                                 y.max.manual=NULL, y.min.manual=NULL, no.title=FALSE, override.title=NULL,
                                 add.chr.to.label=FALSE, alternative.labels=NULL)
{
  allele.effects <- scan.object$allele.effects
  
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
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  allele.effects <- allele.effects[, order.i]
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  
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
  
  plot(0, pch="",
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)), 
       ylim=c(y.min, y.max), 
       xaxt="n", yaxt="n", xlab="", ylab="Additive allele effects", main=this.title,
       frame.plot=FALSE, type="l", cex=0.5, lwd=1.5, col=main.colors[1])
  axis(side=2, at=y.min:y.max, las=2)
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  
  ## Confint for first chr
  if(!is.null(imp.confint.alpha)){
    for(i in 1:num.founders){
      polygon(x=c(pos[pre.chr==chr.types[1]], rev(pos[pre.chr==chr.types[1]])),
              y=c(allele.effect.confint[1, i,  pre.chr==chr.types[1]], rev(allele.effect.confint[2, i, pre.chr==chr.types[1]])),
              border=NA, col=scales::alpha(main.colors[i], transparency))
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
        polygon(x=c(min(this.pos, na.rm=TRUE), min(this.pos, na.rm=TRUE):max(this.pos, na.rm=TRUE), max(this.pos, na.rm=TRUE)), 
                y=c(y.min, rep(y.max, length(min(this.pos, na.rm=TRUE):max(this.pos, na.rm=TRUE))), y.min), border=NA, col="gray88")
      }
      
      ## Confint for later chr
      if(!is.null(imp.confint.alpha)){
        for(j in 1:num.founders){
          polygon(x=c(pos[pre.chr==chr.types[i]] + shift, rev(pos[pre.chr==chr.types[i]] + shift)),
                  y=c(allele.effect.confint[1, j,  pre.chr==chr.types[i]], rev(allele.effect.confint[2, j, pre.chr==chr.types[i]])),
                  border=NA, col=scales::alpha(main.colors[j], transparency))
        }
      }
      ## Means for later chr
      for(j in 1:num.founders){
        points(x=pos[pre.chr==chr.types[i]] + shift, y=allele.effects[j,  pre.chr==chr.types[i]],
               type="l", lwd=my.lwd[j], col=main.colors[j])
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)
      #points(this.pos, allele.effects[1, pre.chr==chr.types[i]], type="l", lwd=1.5, col=main.colors[1])
      shift <- shift + max.pos[i]
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
           fill=main.colors, border="black", 
           col=main.colors, bty=my.bty, cex=my.legend.cex)
  }
}