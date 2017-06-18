
# Ronnegard & Valdar 2011
#' @export
calc.SIM.cache <- function(cache.dir, use.dosages=FALSE, use.subjects=NULL, model="additive", p.null=0.125, scale="Mb", just.these.loci=NULL){
  h <- DiploprobReader$new(cache.dir)
  
  loci <- h$getLoci()
  if(!is.null(just.these.loci)){ loci <- just.these.loci }
  subjects <- h$getSubjects()

  subjects.i <- rep(TRUE, length(subjects))
  if(!is.null(use.subjects)){
    subjects.i <- subjects %in% use.subjects
    subjects <- subjects[subjects.i]
  }
  
  SIM.matrix <- matrix(0, nrow=length(subjects), ncol=length(loci))
  for(i in 1:length(loci)) {
    probs <- h$getLocusMatrix(locus=loci[i], model=model)[subjects.i,]
    if(use.dosages){
      probs <- probs/2
    }
    SIM.vector <- rowSums(get.SIM(p.post=probs, p.null=p.null))
    SIM.matrix[,i] <- SIM.vector
  }
  rownames(SIM.matrix) <- h$getSubjects()[subjects.i]
  colnames(SIM.matrix) <- loci
  return(list(SIM.matrix=SIM.matrix, 
              pos=h$getMarkerLocation(markers=loci, scale=scale),
              chr=h$getChromOfLocus(loci)))
}

get.SIM <- function(p.post, p.null=NULL){
  if(is.null(p.null)){ p.null <- 1/length(p.post) }
  p.post[p.post == 0] <- .Machine$double.eps
  SIM.ind <- p.post*log(p.post/p.null)
  return(SIM.ind/log(1/p.null))
}

#' @export
plot_SIM <- function(SIM.object, chr=1, main=""){
  use.based.chr <- SIM.object$chr == chr
  pos <- SIM.object$pos[use.based.chr]
  SIM <- SIM.object$SIM.matrix[,use.based.chr]
  o <- order(pos)
  plot(x=pos[o], y=SIM[1,o], type="l", col="grey", ylim=c(0, 1), ylab="Selective Information Content (SIC)", xlab=paste0("Chr", chr, " position (Mb)"), main=main, las=1)
  for (i in 2:nrow(SIM))
  {
    lines(x=pos[o], y=SIM.object$SIM[i,o], col="grey")
  }
  lines(x=pos[o], y=colMeans(SIM)[o], col="black", lwd=3)
}

calc.dist.from.balance <- function(genomecache, model){
  require(DiploprobReader)
  h <- DiploprobReader$new(genomecache)
  loci <- h$getLoci()
  dist.vec <- rep(NA, length(loci))
  for(i in 1:length(loci)){
    prob.matrix <- h$getLocusMatrix(loci[i], model=model)
    observed.balance <- colSums(prob.matrix)
    even.balance <- rep(sum(observed.balance)/length(observed.balance), length(observed.balance))
    dist.vec[i] <- as.numeric(dist(rbind(observed.balance, even.balance)))
  }
  dist.vec <- dist.vec/max(dist.vec)
  names(dist.vec) <- loci
  return(dist.vec)
}

#' @export
calc.dist.cache <- function(cache.dir, use.dosages=FALSE, use.subjects=NULL, interval=1) 
{
  h <- DiploprobReader$new(cache.dir)
  
  loci <- h$getLoci()
  o <- order(h$getMarkerLocation(loci, scale="Mb"))
  loci <- loci[o]
  subjects <- h$getSubjects()
  
  subjects.i <- rep(TRUE, length(subjects))
  if (!is.null(use.subjects))
  {
    subjects.i <- subjects %in% use.subjects
    subjects <- subjects[subjects.i]
  }
  
  num.intervals <- floor((length(loci)-1)/interval)
  dist.matrix <- matrix(0, nrow=length(subjects), ncol=num.intervals)
  for (i in 1:length(subjects)) 
  {
    marker.matrix <- matrix(0, nrow=ceiling(length(loci)/interval), ncol=8)
    marker.count <- 1
    for (j in seq(1, length(loci), by=interval)) 
    {
      marker.matrix[marker.count,] <- h$getLocusMatrix(locus=loci[j], model="additive")[subjects.i,][i,]
      marker.count <- marker.count + 1
    }
    if (use.dosages) 
    {
      marker.matrix <- marker.matrix/2
    }
    dist.ind <- diag(as.matrix(dist(marker.matrix))[-1,])
    dist.matrix[i,] <- dist.ind
  }
  rownames(dist.matrix) <- h$getSubjects()[subjects.i]
  colnames(dist.matrix) <- loci[seq(1, length(loci), by=interval)][-(ncol(dist.matrix)+1)]
  return(list(dist.matrix=dist.matrix, 
              Mb=h$getMarkerLocation(markers=loci[seq(1, length(loci), by=interval)][-(ncol(dist.matrix)+1)], scale="Mb"),
              interval.spacing=interval))
}

#' @export
plot_dist <- function(dist.object, comparison.object=NULL, main="", random.sample=20, seed=10){
  o <- order(dist.object$Mb)
  plot(x=dist.object$Mb[o], y=dist.object$dist.matrix[1,o], type="p", pch=20, col="grey", 
       ylim=c(0, max(dist.object$dist.matrix, comparison.object$dist.matrix)), 
       ylab="Euclidian Distance", xlab="Chr1 position (Mb)", 
       main=paste0(main, " - Euclidian distance between neighboring markers (", dist.object$interval, "-marker intervals)"), las=1)
  for (i in 2:nrow(dist.object$dist.matrix)){
    lines(x=dist.object$Mb[o], y=dist.object$dist.matrix[i,o], type="p", pch=20, col="grey")
  }
  set.seed(seed)
  rs <- sample(x=1:nrow(dist.object$dist.matrix), size=random.sample)
  for (i in 1:random.sample){
    lines(x=dist.object$Mb[o], y=dist.object$dist.matrix[random.sample[i],o], type="l", pch=20, col="red")
  }
  lines(x=dist.object$Mb[o], y=colMeans(dist.object$dist.matrix)[o], col="black", lwd=3)
  abline(h=0, lty=2, col="blue")
  abline(h=dist(rbind(c(1,0,0,0,0,0,0,0), c(0,1,0,0,0,0,0,0))), lty=2, col="blue")
}
