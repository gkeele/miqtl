plot.ci <- function(midvals, narrow.intervals, wide.intervals,
                    names=1:length(midvals),
                    add=FALSE,
                    xlab="Estimate", xlab.line=2.5, x.cex=1, x.padj=1, x.labels=TRUE, xlab.cex=1,
                    ylab="", y.cex=1,
                    yaxis=TRUE,
                    name.margin=6,
                    name.line=2,
                    pch.midvals=19, point.cex=1,
                    col=rep("black", length(midvals)),
                    col.midvals=col,
                    shift=0,
                    type="p",
                    use.this.lim=NULL,
                    main="",
                    ...)
{
  nvals <- length(midvals)
  col.midvals <- rep(col.midvals, length.out=nvals)
  y.pos <- (1:nvals)-0.5
  if(!add){
    if(is.null(use.this.lim)){
      lim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
      lim <- c(-1,1) * diff(lim)*0.1 + lim
    }
    if(!is.null(use.this.lim)){
      lim <- use.this.lim
    }
    
    plot(lim, c(0,nvals+0.5), type="n", axes=FALSE, ylab=ylab, xlab="", main=main, ...)
    mtext(text=xlab, side=1, line=xlab.line, cex=xlab.cex)
    axis(1, cex.axis=x.cex, padj=x.padj, labels=x.labels)
    if(yaxis){
      axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line,
           cex.axis=y.cex)
    }
  }
  if("p"==type){
    for(i in 1:nvals){
      pos <- nvals-i + 0.5 + shift
      lines(wide.intervals[i,], rep(pos,2), col=col[i])
      lines(narrow.intervals[i,], rep(pos,2), lwd=3, col=col[i])
      points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i], cex=point.cex)
    }
  }
  invisible(rev(y.pos))
}

#' Plot the regression coefficients as allele effects for a single locus as a caterpillar plot from a scan object
#'
#' This function takes a scan object from scan.h2lmm() and plots out allele effect estimates
#' based on regression coefficients for a single locus. If multiple imputation was used, these correspond to 
#' confidence intervals on the mean effect over imputations.
#' 
#' @param scan.object A scan object from scan.h2lmm(). Can work for ROP or multiple imputation.
#' @param locus The locus for which effect estimates are need. The locus must be included in the scan.
#' @param main DEFAULT: "". The title to included in the plot.
#' @param names DEFAULT: c("ACI", "BN", "BUF", "F344", "M520", "MR", "WKY", "WN"). These are the strains used in the HS rats.
#' @param col DEFAULT: c(rgb(240, 240, 0, maxColorValue=255), rgb(128, 128, 128, maxColorValue=255), rgb(240, 128, 128, maxColorValue=255), 
#' rgb(16, 16, 240, maxColorValue=255), rgb(0, 160, 240, maxColorValue =255),rgb(0, 160, 0, maxColorValue =255),rgb(240, 0, 0, maxColorValue = 255), 
#' rgb(144, 0, 224, maxColorValue = 255)). These are the established Collaborative Cross colors.
#'
#' @export
#' @examples plot.locus.effect.from.scan()
plot.locus.effect.from.scan <- function(scan.object, locus,
                                        main="",
                                        names=c("ACI", "BN", "BUF", "F344", "M520", "MR", "WKY", "WN"),
                                        col=c(rgb(240, 240, 0, maxColorValue=255), 
                                              rgb(128, 128, 128, maxColorValue=255), 
                                              rgb(240, 128, 128, maxColorValue=255), 
                                              rgb(16, 16, 240, maxColorValue=255), 
                                              rgb(0, 160, 240, maxColorValue =255),
                                              rgb(0, 160, 0, maxColorValue =255),
                                              rgb(240, 0, 0, maxColorValue = 255), 
                                              rgb(144, 0, 224, maxColorValue = 255))){
  allele.effects <- scan.object$allele.effects
  effect.matrix <- allele.effects[,locus,]
  
  if(length(col) == 1){ col <- rep(col, nrow(effect.matrix)) }
  
  effects.95 <- t(apply(effect.matrix, 1, function(x) miqtl::ci.mean(x, alpha=0.05)))
  effects.75 <- t(apply(effect.matrix, 1, function(x) miqtl::ci.mean(x, alpha=0.25)))
  mean.effects <- apply(effect.matrix, 1, function(x) mean(x, alpha=0.05, na.rm=TRUE))

  plot.ci(midvals=mean.effects, narrow.intervals=effects.75, wide.intervals=effects.95,
          names=names, col=col)   
  abline(v=0, lty=2)
  title(main)
}