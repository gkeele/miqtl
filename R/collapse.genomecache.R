#' Reduce genome caches by averaging together loci that are very similar up to some set
#' tolerance level of a specified criterion
#'
#' This function takes an inpute genome cache directory, and produces a smaller version
#' with similar information. It is particularly useful for large, dense genome caches, with
#' redundant loci.
#' 
#' @param original.cache The path of the genome cache to be reduced.
#' @param new.cache The path of the new genome cache to be created.
#' @param subjects DEFAULT: NULL. Allows for the specification of a reduced set of individuals from
#' the original genome cache to be included in the new one. The default includes all individuals.
#' @param criterion DEFAULT: "l2.norm". Option to specify criterion for collapsing loci founder 
#' probabilities/dosages. "l2.norm" means that max l2 norm or Euclidean distance is used, which means
#' changes across all categories count. "max.category" means only changes in the max category (founder
#' dosage or diplotype) are used.
#' @param model DEFAULT: "additive". If "additive" is specified, criteria are based on dosages. If 
#' "full", probabilities of diplotypes are used.
#' @param proportion.tol DEFAULT: 0.01. If the maximum criterion value at a pair of loci in a data set 
#' exceeds this value, the loci are not averaged. When all criterion values are below it, the two loci
#' get averaged. The scale of the parameter is in terms of a proportion, regardless of criteria.
#' @export collapse.genomecache
#' @examples collapse.genomecache()
collapse.genomecache <- function(original.cache, 
                                 new.cache, subjects=NULL,
                                 criterion=c("l2.norm", "max.category"),
                                 model=c("additive", "full"),
                                 proportion.tol=0.01){
  criterion <- criterion[1]
  model <- model[1]
  
  h <- DiploprobReader$new(original.cache)
  
  # Chromosomes in original cache
  chr <- list.dirs(paste0(original.cache, "/additive/"), full.names=FALSE, recursive=FALSE)
  chr <- chr[grepl(x=chr, pattern="chr")]
  chr <- gsub(x=chr, pattern="chr", replacement="")

  full.to.add.matrix <- straineff.mapping.matrix()
  
  strains <- h$getFounders()
  if(is.null(subjects)){ subjects <- h$getSubjects() }

  for(i in 1:length(chr)){
    dir.create(paste0(new.cache, "/additive", "/chr", chr[i], "/data"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(new.cache, "/full", "/chr", chr[i], "/data"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(new.cache, "/genotype", "/chr", chr[i]), recursive=TRUE, showWarnings=FALSE)
    
    ## Grab loci from this chr
    load(paste0(original.cache, "/additive", "/chr", chr[i], "/markers.RData"))
    these.loci <- markers
    these.map <- h$getLocusStart(loci=these.loci, scale="cM")

    # Re-ordering based on genetic position
    this.order <- order(these.map)
    these.loci <- these.loci[this.order]
    
    total.loci <- length(these.loci)
    left.marker <- these.loci[1]
    reduced.loci <- reduced.map <- reduced.pos <- NULL
    bin.matrix <- h$getLocusMatrix(locus=left.marker, model="full", subjects=subjects); bin.loci.count <- 1
    for(j in 2:total.loci){
      OUTPUT=FALSE
      right.marker <- these.loci[j]
      
      locus1 <- h$getLocusMatrix(locus=left.marker, model=model, subjects=subjects)
      locus2 <- h$getLocusMatrix(locus=right.marker, model=model, subjects=subjects)
      
      ## Check to see if X changes between markers
      if(criterion == "l2.norm"){
        test <- check.l2.norm(locus1.matrix=locus1, locus2.matrix=locus2, proportion.tol=proportion.tol, model=model)
      }
      else if(criterion == "max.category"){
        test <- check.max.category(locus1.matrix=locus1, locus2.matrix=locus2, proportion.tol=proportion.tol)
      }

      ## Extend bin
      if(test){
        bin.matrix <- bin.matrix + h$getLocusMatrix(locus=right.marker, model="full", subjects=subjects)
        bin.loci.count <- bin.loci.count + 1
        ## CASE: Muli-locus bin at end of chromosome
        if(j == total.loci){
          bin.marker <- paste0(left.marker, "_to_", right.marker) # Bin names
          bin.matrix <- bin.matrix/bin.loci.count # Averaging markers
          OUTPUT=TRUE
        }
      }
      
      ## End bin and start new one
      if(!test){
        ## Bin is a single marker
        if(bin.loci.count == 1){
          bin.marker <- left.marker
        }
        else{
          bin.marker <- paste0(left.marker, "_to_", these.loci[j-1]) ## Bin names
          bin.matrix <- bin.matrix/bin.loci.count ## Averaging markers
        }
        OUTPUT=TRUE
      }
      
      #browser()
      ## Adding bin to output
      if(OUTPUT){
        reduced.map <- c(reduced.map, h$getLocusStart(loci=left.marker, scale="cM"))
        reduced.pos <- c(reduced.pos, h$getLocusStart(loci=left.marker, scale="Mb"))
        reduced.loci <- c(reduced.loci, bin.marker)
        
        # Additive dosages
        bin.dosages <- bin.matrix %*% full.to.add.matrix
        colnames(bin.dosages) <- strains
        rownames(bin.dosages) <- rownames(bin.matrix)
        assign(bin.marker, bin.dosages)
        add.fn <- paste0(new.cache, "/additive/chr", chr[i], "/data/", bin.marker, ".RData")
        save(list=bin.marker, file=add.fn)
        
        # Full probabilities
        assign(bin.marker, bin.matrix)
        full.fn <- paste0(new.cache, "/full/chr", chr[i], "/data/", bin.marker, ".RData")
        save(list=bin.marker, file=full.fn)
        
        # Removing locus matrix so memory doesn't get overloaded
        rm(list=bin.marker)
      }
      
      ## CASE: Last marker is alone in a bin
      if(!test & j == total.loci){
        bin.marker <- right.marker
        bin.matrix <- h$getLocusMatrix(locus=bin.marker, model="full", subjects=subjects)
        
        # Additive dosages
        bin.dosages <- bin.matrix %*% full.to.add.matrix
        colnames(bin.dosages) <- strains
        rownames(bin.dosages) <- rownames(bin.matrix)
        assign(bin.marker, bin.dosages)
        add.fn <- paste0(new.cache, "/additive/chr", chr[i], "/data/", bin.marker, ".RData")
        save(list=bin.marker, file=add.fn)
        
        # Full probabilities
        assign(bin.marker, bin.matrix)
        full.fn <- paste0(new.cache, "/full/chr", chr[i], "/data/", bin.marker, ".RData")
        save(list=bin.marker, file=full.fn)
        
        # Removing locus matrix so memory doesn't get overloaded
        rm(list=bin.marker)
        
        reduced.map <- c(reduced.map, h$getLocusStart(loci=right.marker, scale="cM"))
        reduced.pos <- c(reduced.pos, h$getLocusStart(loci=right.marker, scale="Mb"))
        reduced.loci <- c(reduced.loci, right.marker)
      }
        
      # Reset
      if(!test & j < total.loci){
        left.marker <- right.marker 
        bin.matrix <- h$getLocusMatrix(locus=left.marker, model="full", subjects=subjects) 
        bin.loci.count <- 1
      }
    }
    ## Markers
    markers <- reduced.loci
    save(list="markers", file=paste0(new.cache, "/full/chr", chr[i], "/markers.RData"))
    save(list="markers", file=paste0(new.cache, "/additive/chr", chr[i], "/markers.RData"))
    save(list="markers", file=paste0(new.cache, "/genotype/chr", chr[i], "/markers.RData"))
    
    ## Map
    #map <- unique(reduced.map)
    map <- reduced.map
    save(list="map", file=paste0(new.cache, "/full/chr", chr[i], "/map.RData"))
    save(list="map", file=paste0(new.cache, "/additive/chr", chr[i], "/map.RData"))
    save(list="map", file=paste0(new.cache, "/genotype/chr", chr[i], "/map.RData"))
    
    ## Pos
    #bp <- unique(reduced.pos*1e6)
    bp <- reduced.pos*1e6
    save(list="bp", file=paste0(new.cache, "/full/chr", chr[i], "/bp.RData"))
    save(list="bp", file=paste0(new.cache, "/additive/chr", chr[i], "/bp.RData"))
    save(list="bp", file=paste0(new.cache, "/genotype/chr", chr[i], "/bp.RData"))
    
    ## Chromosome
    chromosome <- rep(chr[i], length(reduced.loci))
    save(list="chromosome", file=paste0(new.cache, "/full/chr", chr[i], "/chromosome.RData"))
    save(list="chromosome", file=paste0(new.cache, "/additive/chr", chr[i], "/chromosome.RData"))
    save(list="chromosome", file=paste0(new.cache, "/genotype/chr", chr[i], "/chromosome.RData"))
    
    ## Strains
    save(list="strains", file=paste0(new.cache, "/full/chr", chr[i], "/strains.RData"))
    save(list="strains", file=paste0(new.cache, "/additive/chr", chr[i], "/strains.RData"))
    save(list="strains", file=paste0(new.cache, "/genotype/chr", chr[i], "/strains.RData"))
    
    ## Subjects
    save(list="subjects", file=paste0(new.cache, "/full/chr", chr[i], "/subjects.RData"))
    save(list="subjects", file=paste0(new.cache, "/additive/chr", chr[i], "/subjects.RData"))
    save(list="subjects", file=paste0(new.cache, "/genotype/chr", chr[i], "/subjects.RData"))
    
    cat(paste0("Finished Chr", chr[i], "\n"))
  }
}

## Returns TRUE if l2norm for all individuals is below some tolerance level
check.l2.norm <- function(locus1.matrix, locus2.matrix, proportion.tol, model){
  dif.matrix <- locus1.matrix - locus2.matrix
  
  max.l2 <- ifelse(model=="additive", sqrt(8), sqrt(2))
  
  l2.norm <- apply(dif.matrix, 1, function(x) sqrt(sum(x^2))/max.l2)
  
  return(ifelse(any(l2.norm > proportion.tol), FALSE, TRUE))
}

## Returns TRUE if difference in max category for all individuals is below some tolerance level
## In practice, a max haplotype dosage or max diplotype
check.max.category <- function(locus1.matrix, locus2.matrix, proportion.tol){
  max.category <- apply(locus1.matrix, 1, function(x) which.max(x))
  
  cat.dif <- sapply(1:length(max.category), 
                    function(x) abs(locus1.matrix[x, max.category[x]] - locus2.matrix[x, max.category[x]]))
  
  return(ifelse(any(cat.dif > proportion.tol), FALSE, TRUE))
}




