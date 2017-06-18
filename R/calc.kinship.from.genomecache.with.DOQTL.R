#' @export
calc.kinship.from.genomecache.with.DOQTL <- function(genomecache, model="additive"){
  h <- DiploprobReader$new(genomecache)
  mapping <- straineff.mapping.matrix()
  
  probs <- array(NA, dim=c(length(h$getSubjects()), length(h$getFounders()), length(h$getLoci())))
  
  loci <- h$getLoci()
  
  for(i in 1:length(loci)){
    diplotypes <- h$getLocusMatrix(loci[i], model="full")
    if(any(diplotypes < 0)){
      diplotypes[diplotypes < 0] <- 0
      diplotypes <- t(apply(diplotypes, 1, function(x) x/sum(x)))
    }
    if(model == "additive"){
      use.probs <- (diplotypes %*% mapping)/2
    }
    if(model == "full"){
      use.probs <- diplotypes
    }
    probs[,,i] <- use.probs
  }
  
  K <- DOQTL::kinship.probs(probs)
  colnames(K) <- rownames(K) <- h$getSubjects()
  return(K)
}