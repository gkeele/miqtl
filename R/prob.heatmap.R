#' Plot founder haplotype dosages/probabilities, ordered by phenotype
#'
#' This function produces a probability heatmap plot, ordered by the phenotype. This gives an idea of what the regression
#' procedure is actually being handed as inputs.
#'
#' @param marker A marker that is contained in the genome cache. In general this should be a marker of interest, such as one
#' beneath a putative QTL peak.
#' @param p.value DEFAULT: NULL. Includes the observed p-value in the plot title.
#' @param genomecache The path to the genome cache that contains founder haplotype information.
#' @param model DEFAULT: "additive". If "additive", dosages are plotted. If "full", probabilities are plotted.
#' @param phenotype The name of the phenotype column in the data set.
#' @param phenotype.data A data.frame object that contains the phenotype information. Should also have a column that matches
#' genomes in the genome cache.
#' @param merge.by DEFAULT: "SUBJECT.NAME". Specifies the columns to merge phenotype and haplotype data.
#' @param founder.labels DEFAULT: "NULL". If NULL, will default to the labels in the genome cache.
#' @param founder.cex DEFAULT: 1. Defines the text size of the founder labels.
#' @param include.ramp DEFAULT: TRUE. If TRUE, spectrum ramp for dosages or probabilities is included.
#' @param ramp.label.cex DEFAULT: 0.7. Specifies the size of the dosage/probability spectrum label.
#' @param ramp.label.line DEFAULT: 0.5. Specifies the number of margin lines that separates the spectrum label from the bottom axis.
#' @param prob.axis.cex DEFAULT: 1. Specifies the size of the dosage/probability tick labels.
#' @param include.marker DEFAULT: TRUE. If TRUE, the marker name is included in the title.
#' @param alternative.phenotype.label DEFAULT: NULL. Allows for an alternative label for the phenotype.
#' @param alternative.marker.label DEFAULT: NULL. Allows for an alternative label for the marker.
#' @export
#' @examples prob.heatmap()
prob.heatmap = function(marker, p.value=NULL, genomecache, model="additive",
                        phenotype, phenotype.data, merge.by="SUBJECT.NAME", 
                        founder.labels=NULL, founder.cex=1,
                        phenotype.lab.cex=1, phenotype.num.cex=1, phenotype.num.padj=NA,
                        phenotype.line=NA, phenotype.num.line=NA,
                        include.ramp=TRUE, ramp.label.cex=0.7, ramp.label.line=0.5, prob.axis.cex=1,
                        include.marker=TRUE, marker.line=1,
                        alternative.phenotype.label=NULL, alternative.marker.label=NULL){
  h <- DiploprobReader$new(genomecache)
  X <- h$getLocusMatrix(locus=marker, model=model)
  subjects <- h$getSubjects()
  rownames(X) <- subjects
  
  if(!include.marker){ marker <- NULL }
  if(!is.null(alternative.marker.label)){ marker <- alternative.marker.label }
  
  prob.heatmap.from.matrix(geno.matrix=X, marker=marker, p.value=p.value, model=model, phenotype=phenotype,
                           phenotype.data, merge.by=merge.by, 
                           founder.labels=founder.labels, founder.cex=founder.cex,
                           phenotype.lab.cex=phenotype.lab.cex, phenotype.num.cex=phenotype.num.cex, phenotype.num.padj=phenotype.num.padj,
                           phenotype.line=phenotype.line, phenotype.num.line=phenotype.num.line,
                           include.ramp=include.ramp, ramp.label.cex=ramp.label.cex, ramp.label.line=ramp.label.line, prob.axis.cex=prob.axis.cex,
                           marker.line=marker.line,
                           alternative.phenotype.label=alternative.phenotype.label)
}

#' @export
prob.heatmap.from.matrix = function(geno.matrix, marker=NULL, marker.line=1,
                                    p.value=NULL, model="additive",
                                    phenotype, phenotype.data,
                                    merge.by="SUBJECT.NAME", 
                                    founder.labels=NULL, founder.cex=1, 
                                    phenotype.lab.cex=1, phenotype.num.cex=1, phenotype.num.padj=NA,
                                    phenotype.line=NA, phenotype.num.line=NA,
                                    include.ramp=TRUE, ramp.label.cex=0.7, ramp.label.line=0.5, prob.axis.cex=1,
                                    alternative.phenotype.label=NULL){
  if(!is.null(p.value)){
    p.value <- round(-log10(p.value), 4)
  }
  
  X <- geno.matrix
  X.data <- data.frame(rownames(X), X)
  names(X.data)[1] <- merge.by
  
  if(is.null(founder.labels)){
    founder.labels <- colnames(X)
  }
  num.col <- length(founder.labels)
  
  # Allow function of phenotype
  phenotype.data <- model.frame(formula(paste(phenotype, "~ 1 +", merge.by)), data=phenotype.data)
  names(phenotype.data)[1] <- "y"
  final.data <- merge(x=phenotype.data, y=X.data, by=merge.by, all=FALSE)
  
  final.data <- final.data[order(final.data$y),] # sort by phenotypic value
  probs <- as.matrix(final.data[,-(1:2)]) # Just keep prob of 8 strains for a certain marker
  probs <- probs[,rev(1:ncol(probs))]
  
  s <- summary(final.data[,"y"])
  s1 <- as.character(round(s[1], 2))
  s2 <- as.character(round(s[2], 2))
  s3 <- as.character(round(s[3], 2))
  s5 <- as.character(round(s[5], 2))
  s6 <- as.character(round(s[6], 2))
  
  # plot heatmap
  ## save original par settings
  cols <- rev(gray(10000:1/10000))
  if(include.ramp){
    op <- par(no.readonly=TRUE)
    oplt <- par()$plt
    par(fig=c(0.1, 0.85, 0.05, 0.95),
        mai=c(0.4, 0.05, 0.7, 0.05))    ##set the margin  
  }
  if(model == "additive"){ z.val <- 2 - probs; z.lim <- c(0, 2) }
  else{ z.val <- 1 - probs; z.lim <- c(0, 1) }
  image(z=z.val, axes=FALSE, col=cols, zlim=z.lim)
  box()
  axis(2, at=seq(0, num.col, 1+1/num.col)/num.col, 
       labels=rev(founder.labels), cex.axis=founder.cex,
       lty=0, srt=90, las=2) # add txt on the strain 
  phenotype <- ifelse(is.null(alternative.phenotype.label), phenotype, alternative.phenotype.label)
  axis(1, at=0.5, labels=parse(text=paste('"-" %<-%', phenotype, '%->% "+"')), tick=FALSE, cex.axis=phenotype.lab.cex, line=phenotype.line)
  axis(3, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(s1, s2, s3, s5, s6), cex.axis=phenotype.num.cex, line=phenotype.num.line, padj=phenotype.num.padj)
  if(!is.null(marker)){
    if(is.null(p.value)){
      this.title <- marker
    }
    else{
      this.title <- bquote(.(paste0(marker, ":")) ~ -log[10]*P ~ .(paste("=", p.value)))
    }
    title(this.title, line=marker.line)
  }
  
  if(include.ramp){
    ramp.label <- ifelse(model == "additive", "Haplotype Dose", "Diplotype Prob")
    par(fig=c(0.89, 0.92, 0.33, 0.66), 
        mai=c(0.1, 0.05, 0.5, 0.05), 
        new=TRUE)
    if(model == "additive"){ image(y=seq(from=0, to=2, length.out=length(cols)), z=matrix(seq(from=0, to=2, length.out=length(cols)), nrow=1), 
                                   zlim=c(0, 2), ylim=c(0, 2), axes=FALSE, col=rev(cols), main="", cex.main=0.77) } #for the legend 
    if(model == "full"){ image(y=seq(from=0, to=1, length.out=length(cols)), z=matrix(seq(from=0, to=1, length.out=length(cols)), nrow=1), 
                               zlim=c(0, 1), ylim=c(0, 1), axes=FALSE, col=rev(cols), main="", cex.main=0.77) }
    box()
    axis(4, las=1, cex.axis=prob.axis.cex)
    mtext(text=ramp.label, side=1, line=ramp.label.line, cex=ramp.label.cex)
    par(op)
  }
  #else{ par(plt <- oplt) }
}


#' @export
prob.image = function(marker.data, marker=NULL, p.value=NULL, 
                      phenotype, phenotype.data, column.labels=NULL){
  if(!is.null(p.value)){ p.value <- round(-log10(p.value), 4) }
  if(is.null(column.labels)){ column.labels <- colnames(marker.data)[-1] }
  if(is.null(marker)){ marker <- "locus" }
  num.col <- length(column.labels)
  if(sum(marker.data[1,-1], na.rm=TRUE) != 1){
    marker.data[,-1] <- marker.data[,-1]/2
  }
  # allow function of phenotype
  phenotype.data <- model.frame(formula(paste0(phenotype, " ~ 1 + SUBJECT.NAME")), data=phenotype.data)
  names(phenotype.data)[1] <- "y"
  final.data <- merge(x=phenotype.data, y=marker.data, by="SUBJECT.NAME", all=FALSE)
  
  final.data <- final.data[order(final.data$y),] # sort by phenotypic value
  probs <- as.matrix(final.data[,-(1:2)]) # Just keep prob of 8 strains for a certain marker
  
  s <- summary(final.data[,"y"])
  s1 <- as.character(round(s[1], 2))
  s2 <- as.character(round(s[2], 2))
  s3 <- as.character(round(s[3], 2))
  s5 <- as.character(round(s[5], 2))
  s6 <- as.character(round(s[6], 2))
  
  # plot image
  ## save original par settings
  oplt <- par()$plt
  cols <- rev(gray(10000:1/10000))
  par(plt=c(0.1,.75,0.1,.8))    ##set the margin  
  image(z=1-probs, axes=F, col=cols, zlim=c(0, 1))
  box()
  axis(2, at=seq(0, num.col, 1+1/num.col)/num.col, 
       labels=column.labels, 
       lty=0, srt=90, las=2) # add txt on the strain   
  axis(1, at=0.5, labels=phenotype, tick=FALSE)
  axis(3, at=c(0,0.25,0.5,0.75,1), labels=c(s1,s2,s3,s5,s6))
  this.title <- ifelse(is.null(p.value), marker, paste0(marker, ": -log10P=", p.value))
  title(this.title, line=2.5)  
  par(plt <- oplt)
}
