convert.qtl2.to.HAPPY <- function(qtl2.object, map,
                                  map.locus_name.colname="SNP_ID", map.chr.colname="Chr", map.physical_dist.colname="Mb_NCBI38", map.genetic_dist.colname="cM",
                                  HAPPY.output.path, remove.chr.from.chr=FALSE,
                                  physical_dist.is.Mb=TRUE,
                                  allele.labels=LETTERS[1:8],
                                  chr=c(1:19, "X"),
                                  diplotype.order=c("qtl2", "DOQTL", "CC")){
  
  full.array.list <- qtl2.object$probs
  total.map <- map
  
  samples <- dimnames(full.array.list[[1]])[[1]]
  diplotypes <- dimnames(full.array.list[[1]])[[2]]
  diplotype.order <- diplotype.order[1]
  
  chr <- names(full.array.list)[names(full.array.list) %in% chr]
  
  full.to.add.matrix <- straineff.mapping.matrix()
  ## Reorder diplotype columns
  if(diplotype.order == "qtl2"){
    dip.order <- c(1, 3, 6, 10, 15, 21, 28, 36,
                   2, 
                   4, 5,
                   7, 8, 9,
                   11, 12, 13, 14,
                   16, 17, 18, 19, 20,
                   22, 23, 24, 25, 26, 27,
                   29, 30, 31, 32, 33, 34, 35)
  }
  else if(diplotype.order == "DOQTL"){
    dip.order <- c(1,9,16,22,27,31,34,36,2,3,10,4,11,
                   17,5,12,18,23,6,13,19,24,28,7,14,
                   20,25,29,32,8,15,21,26,30,33,35)
  }
  else if(diplotype.order == "CC"){
    dip.order <- c(1,2,3,4,5,6,7,8,
                   9,
                   10,16,
                   11,17,22,
                   12,18,23,27,
                   13,19,24,28,31,
                   14,20,25,29,32,34,
                   15,21,26,30,33,35,36)
  }
  
  ## Making necessary directories
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/additive/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/full/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/genotype/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  
  subjects <- samples
  strains <- allele.labels
  for(i in 1:length(chr)){
    loci <- dimnames(full.array.list[[chr[i]]])[[3]]
    
    ## Reducing map to just those also in array
    chr.total.map <- total.map[total.map[,map.locus_name.colname] %in% loci,]
    ## Reducing loci to just those also in map
    loci <- loci[loci %in% chr.total.map[,map.locus_name.colname]]

    for(j in 1:length(loci)){
      chr.locus <- as.character(chr.total.map[chr.total.map[,map.locus_name.colname] == loci[j], map.chr.colname])
      
      locus.matrix <- full.array.list[[chr[i]]][,dip.order,loci[j]]
      
      ## Handling qtl2's approach to X
      if(chr.locus == "X"){
        Y.matrix <- full.array.list[[chr[i]]][,-(1:36),loci[j]]
        
        locus.matrix[,"AA"] <- locus.matrix[,"AA"] + Y.matrix[,"AY"]
        locus.matrix[,"BB"] <- locus.matrix[,"BB"] + Y.matrix[,"BY"]
        locus.matrix[,"CC"] <- locus.matrix[,"CC"] + Y.matrix[,"CY"]
        locus.matrix[,"DD"] <- locus.matrix[,"DD"] + Y.matrix[,"DY"]
        locus.matrix[,"EE"] <- locus.matrix[,"EE"] + Y.matrix[,"EY"]
        locus.matrix[,"FF"] <- locus.matrix[,"FF"] + Y.matrix[,"FY"]
        locus.matrix[,"GG"] <- locus.matrix[,"GG"] + Y.matrix[,"GY"]
        locus.matrix[,"HH"] <- locus.matrix[,"HH"] + Y.matrix[,"HY"]
        
        locus.matrix <- locus.matrix[,1:36]
      }
      
      var_name <- loci[j]
      assign(var_name, locus.matrix)
      full.fn <- paste0(HAPPY.output.path, "/full/chr", chr.locus, "/data/", var_name, ".RData")
      save(list=var_name, file=full.fn)

      dosage.matrix <- locus.matrix %*% full.to.add.matrix
      colnames(dosage.matrix) <- allele.labels
      rownames(dosage.matrix) <- rownames(locus.matrix)
      assign(var_name, dosage.matrix)
      add.fn <- paste0(HAPPY.output.path, "/additive/chr", chr.locus, "/data/", var_name, ".RData")
      save(list=var_name, file=add.fn)
      rm(list=var_name)
    }
    ## chr
    chromosome <- rep(as.character(chr[i]), length(loci))
    save(chromosome, file = paste0(HAPPY.output.path, '/additive/chr', chr[i], '/chromosome.RData'))
    save(chromosome, file = paste0(HAPPY.output.path, '/full/chr', chr[i], '/chromosome.RData'))
    save(chromosome, file = paste0(HAPPY.output.path, '/genotype/chr', chr[i], '/chromosome.RData'))
    ## marker
    markers <- as.character(loci)
    save(markers, file = paste0(HAPPY.output.path, '/additive/chr', chr[i], '/markers.RData'))
    save(markers, file = paste0(HAPPY.output.path, '/full/chr', chr[i], '/markers.RData'))
    save(markers, file = paste0(HAPPY.output.path, '/genotype/chr', chr[i], '/markers.RData'))
    ## bp
    bp <- chr.total.map[, map.physical_dist.colname]
    if(physical_dist.is.Mb){ bp <- bp*1000000 }
    save(bp, file = paste0(HAPPY.output.path, '/additive/chr', chr[i], '/bp.RData'))
    save(bp, file = paste0(HAPPY.output.path, '/full/chr', chr[i], '/bp.RData'))
    save(bp, file = paste0(HAPPY.output.path, '/genotype/chr', chr[i], '/bp.RData'))
    ## map
    map <- chr.total.map[, map.genetic_dist.colname]
    save(map, file = paste0(HAPPY.output.path, '/additive/chr', chr[i], '/map.RData'))
    save(map, file = paste0(HAPPY.output.path, '/full/chr', chr[i], '/map.RData'))
    save(map, file = paste0(HAPPY.output.path, '/genotype/chr', chr[i], '/map.RData'))
    ## subjects
    save(subjects, file = paste0(HAPPY.output.path, '/additive/chr', chr[i], '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/full/chr', chr[i], '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/genotype/chr', chr[i], '/subjects.RData'))
    ## strains
    save(strains, file = paste0(HAPPY.output.path, '/additive/chr', chr[i], '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/full/chr', chr[i], '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/genotype/chr', chr[i], '/strains.RData'))
  }
}
  