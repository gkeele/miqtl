#' Takes DO-QTL founder haplotype reconstruction output files and re-formats into a HAPPY genome cache
#'
#' This function produces a HAPPY-format genome cache from DO-QTL founder haplotype reconstruction output files.
#' The main output files of importance are the individual-level files with naming scheme [sample].genotype.probs.Rdata.
#'
#' @param DOQTL.recon.output.path The path to the directory containing DO-QTL founder haplotype output files. 
#' @param map The map file (which contains important information on the loci) loaded into R.
#' @param physical_dist.is.Mb DEFAULT: TRUE IF true, Mb will be converted to bp within the function.
#' @param map.locus_name.colname DEFAULT: "SNP_ID". The column name in the map data that corresponds to locus/marker names.
#' @param map.chr.colname DEFAULT: "Chr". The column name in the map data that corresponds to chromosome.
#' @param map.physical_dist.colname DEFAULT: "Mb_NCBI38". The column name in the map data that corresponds to physical position.
#' @param map.genetic_dist.colname DEFAULT: "cM". The column name in the map data that corresponds to genetic position.
#' @param HAPPY.output.path The path to a directory that will be created as the HAPPY-format genome cache.
#' @param allele.labels DEFAULT: NULL. Allows for specification of founder labels different from what is in the DO-QTL
#' output. The DEFAULT of NULL leads to using the labels from the DO-QTL output.
#' @param chr DEFAULT: c(1:19, "X"). Allows for specification of the chromosomes. DEFAULT is all the chromosomes from the mouse.
#' @export
#' @import data.table
#' @examples convert.DOQTL.to.HAPPY()
convert.DOQTL.to.HAPPY <- function(DOQTL.recon.output.path,
                                   map, physical_dist.is.Mb=TRUE, 
                                   map.locus_name.colname="SNP_ID", map.chr.colname="Chr", map.physical_dist.colname="Mb_NCBI38", map.genetic_dist.colname="cM",
                                   HAPPY.output.path,
                                   allele.labels=NULL,
                                   chr=c(1:19, "X"),
                                   in.log.scale=FALSE){
  
  #require(data.table)
  #----------------------------------
  # founder probs from DO-QTL
  #----------------------------------
  load(paste(DOQTL.recon.output.path, "founder.probs.Rdata", sep="/"))
  
  samples <- dimnames(model.probs)[[1]]
  simple.alleles <- dimnames(model.probs)[[2]]
  rm(model.probs)
  
  
  #----------------------------------
  # Putting together strain labels
  #----------------------------------
  if(is.null(allele.labels)){
    allele.labels <- simple.alleles
  }

  full.to.dosages <- straineff.mapping.matrix()
  
  diplotype.labels <- c(paste(allele.labels, allele.labels, sep="."), apply(full.to.dosages[-(1:8),], 1, 
                                                                            function(x) paste(allele.labels[sort(which(x==1), decreasing=FALSE)], collapse=".")))
  simple.het.labels <- c(paste(simple.alleles, simple.alleles, sep=""), apply(full.to.dosages[-(1:8),], 1, 
                                                                              function(x) paste(simple.alleles[sort(which(x==1), decreasing=FALSE)], collapse="")))
  
  
  #-------------------------------
  # Marker info
  #-------------------------------
  #map <- read.table(map.path, header=TRUE, as.is=TRUE)
  
  #-------------------------------
  # Functions to output marker files
  #-------------------------------
  output.marker.file <- function(chr, marker,
                aa, bb, cc, dd, ee, ff, gg, hh,
                ba, ca, cb, da, db, dc, ea, eb, ec, ed,
                fa, fb, fc, fd, fe, ga, gb, gc, gd, ge, gf,
                ha, hb, hc, hd, he, hf, hg,
                allele.labels,
                diplotype.labels,
                full.to.dosages.matrix){
    var_name <- marker[1]
    assign(var_name, matrix(data=c(aa, bb, cc, dd, ee, ff, gg, hh,
                                   ba, ca, cb, da, db, dc, ea, eb, ec, ed,
                                   fa, fb, fc, fd, fe, ga, gb, gc, gd, ge, gf,
                                   ha, hb, hc, hd, he, hf, hg),
                            ncol=length(diplotype.labels),
                            dimnames=list(NULL, diplotype.labels)))

    temp <- get(var_name)
    colnames(temp) <- diplotype.labels
    temp.add <- temp %*% full.to.dosages.matrix
    colnames(temp.add) <- allele.labels

    dir.create(paste0(HAPPY.output.path, '/full/chr', chr[1], '/data/'),
               showWarnings=FALSE, recursive=TRUE)
    fn <- paste0(HAPPY.output.path, '/full/chr', chr[1], '/data/', var_name, '.RData')
    save(list=var_name, file=fn)

    assign(var_name, temp.add)
    dir.create(paste0(HAPPY.output.path, '/additive/chr', chr[1], '/data/'),
               showWarnings=FALSE, recursive=TRUE)
    fn <- paste0(HAPPY.output.path, '/additive/chr', chr[1], '/data/', var_name, '.RData')
    save(list = var_name, file = fn)
  }
  # Export file of marker names for each chr
  export_marker_name_file <- function(chr, marker) {
    markers <- as.character(marker)
    save(markers, file=paste0(HAPPY.output.path, '/additive/chr', chr[1], '/markers.RData'))
    save(markers, file=paste0(HAPPY.output.path, '/full/chr', chr[1], '/markers.RData'))
    dir.create(paste0(HAPPY.output.path, '/genotype/chr', chr[1]), showWarnings=FALSE, recursive=TRUE)
    save(markers, file=paste0(HAPPY.output.path, '/genotype/chr', chr[1], '/markers.RData'))
  }
  # Export file of marker bp positions for each chr
  export_marker_position_file <- function(chr, pos) {
    bp <- as.character(pos)
    save(bp, file=paste0(HAPPY.output.path, '/additive/chr', chr[1], '/bp.RData'))
    save(bp, file=paste0(HAPPY.output.path, '/full/chr', chr[1], '/bp.RData'))
    save(bp, file=paste0(HAPPY.output.path, '/genotype/chr', chr[1], '/bp.RData'))
  }
  # Export file of chromosome
  export_marker_chromosome_file <- function(chr, marker) {
    chromosome <- rep(as.character(chr), length(marker))
    save(chromosome, file=paste0(HAPPY.output.path, '/additive/chr',chr[1], '/chromosome.RData'))
    save(chromosome, file=paste0(HAPPY.output.path, '/full/chr',chr[1], '/chromosome.RData'))
    save(chromosome, file=paste0(HAPPY.output.path, '/genotype/chr',chr[1], '/chromosome.RData'))
  }
  # Export file of map distance (cM)
  export_marker_map_distance_file <- function(chr, pos) {
    map <- as.character(pos)
    save(map, file=paste0(HAPPY.output.path, '/additive/chr', chr[1], '/map.RData'))
    save(map, file=paste0(HAPPY.output.path, '/full/chr', chr[1], '/map.RData'))
    save(map, file=paste0(HAPPY.output.path, '/genotype/chr', chr[1], '/map.RData'))
  }
  
  #-----------------------------
  # Rename map data columns
  #-----------------------------
  colnames(map)[colnames(map) == map.locus_name.colname] <- "marker"
  colnames(map)[colnames(map) == map.chr.colname] <- "chr"
  colnames(map)[colnames(map) == map.physical_dist.colname] <- "bp"
  colnames(map)[colnames(map) == map.genetic_dist.colname] <- "pos"
  map <- map[,c("marker", "chr", "bp", "pos")]
  
  #-----------------------------
  # Combining data of individuals
  #-----------------------------
  for(i in 1:length(chr)){
    for(j in 1:length(samples)){
      cat(paste("Loading DOQTL output for individual", j, "for chr", chr[i]), "\n")
      load(paste0(DOQTL.recon.output.path, "/", samples[j], ".genotype.probs.Rdata"))
      if(in.log.scale){
        prsmth <- exp(prsmth)
      }
      marker <- rownames(prsmth)
      subject <- rep(samples[j], nrow(prsmth))
      one.sample.data <- data.frame(marker, subject,  prsmth)
      combined.data <- merge(x=map, y=one.sample.data, by.x="marker", by.y="marker")[,c(1:5,c(1,9,16,22,27,31,34,36,2,3,10,4,11,
                                                                                              17,5,12,18,23,6,13,19,24,28,7,14,
                                                                                              20,25,29,32,8,15,21,26,30,33,35)+5)]
      if(physical_dist.is.Mb){ combined.data$bp <- combined.data$bp*1000000 }
      combined.data <- combined.data[combined.data$chr == chr[i],]
    
      if(!exists('all.subjects')){
        all.subjects <- combined.data
      } 
      else{ 
        all.subjects <- data.table::rbindlist(list(all.subjects, combined.data))
      }
    }
    #--------------------------------------------------------------------------
    # make each marker_name.Rdata
    # Subject order in marker_name.Rdata should match with SUBJECT.NAME in pheno
    #---------------------------------------------------------------------------
    var.names <- c(names(all.subjects)[1:5], simple.het.labels)
    data.table::setnames(all.subjects, names(all.subjects), var.names)
    data.table::setkey(all.subjects, NULL)
    data.table::setkey(all.subjects, chr, bp, marker, subject)

    all.subjects[, output.marker.file(chr, marker,
                                      AA, BB, CC, DD,
                                      EE, FF, GG, HH,
                                      AB, AC, BC, AD,
                                      BD, CD, AE, BE,
                                      CE, DE, AF, BF,
                                      CF, DF, EF, AG,
                                      BG, CG, DG, EG,
                                      FG, AH, BH, CH,
                                      DH, EH, FH, GH,
                                      allele.labels, diplotype.labels, full.to.dosages), by="marker"]
    
    #-------------------------------------
    # make other necessary files
    #-------------------------------------
    one.subj = samples[1]
    markers.one.subj <- all.subjects[grepl(one.subj, subject), ]  
    markers.one.subj[, export_marker_name_file(chr, marker), by="chr"]
    markers.one.subj[, export_marker_position_file(chr, bp), by="chr"]
    markers.one.subj[, export_marker_chromosome_file(chr, marker), by="chr"]
    markers.one.subj[, export_marker_map_distance_file(chr, pos), by="chr"]
    
    rm(all.subjects)
  }
  
  # export file of mice names and strains
  subjects <- samples
  strains <- allele.labels
  for(this.chr in chr){
    save(subjects, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/subjects.RData'))
    
    save(strains, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/strains.RData'))
  }
}

#' @export
convert.additive.DOQTL.array.to.HAPPY <- function(DOQTL.array, map,
                                                  map.locus_name.colname="SNP_ID", map.chr.colname="Chr", map.physical_dist.colname="Mb_NCBI38", map.genetic_dist.colname="cM",
                                                  HAPPY.output.path,
                                                  physical_dist.is.Mb=TRUE,
                                                  allele.labels=NULL, convert.to.dosage=TRUE,
                                                  chr=c(1:19, "X")){
  
  samples <- dimnames(DOQTL.array)[[1]]
  simple.alleles <- dimnames(DOQTL.array)[[2]]
  loci <- dimnames(DOQTL.array)[[3]]
  
  
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/additive/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/full/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/genotype/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  
  #----------------------------------
  # Putting together strain labels
  #----------------------------------
  if(is.null(allele.labels)){
    allele.labels <- simple.alleles
  }
  
  #-------------------------------
  # Marker info
  #-------------------------------
  total.map <- map
  
  for(i in 1:length(loci)){
    chr.locus <- as.character(total.map[total.map$SNP_ID == loci[i], map.chr.colname])
    
    if(convert.to.dosage){ locus.matrix <- DOQTL.array[,,i]*2 }
    else{ locus.matrix <- DOQTL.array[,,i] }
    
    var_name <- loci[i]
    assign(var_name, locus.matrix)
    fn <- paste0(HAPPY.output.path, "/additive/chr", chr.locus, "/data/", var_name, ".RData")
    save(list=var_name, file=fn)
  }
  
  subjects <- samples
  strains <- allele.labels
  for(this.chr in chr){
    ## chr
    map.this.chr <- total.map[as.character(total.map[, map.chr.colname]) == this.chr,]
    chromosome <- as.character(map.this.chr[,map.chr.colname])
    save(chromosome, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/chromosome.RData'))
    save(chromosome, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/chromosome.RData'))
    save(chromosome, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/chromosome.RData'))
    
    ## marker
    markers <- as.character(map.this.chr[, map.locus_name.colname])
    save(markers, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/markers.RData'))
    save(markers, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/markers.RData'))
    save(markers, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/markers.RData'))
    
    ## bp
    bp <- map.this.chr[, map.physical_dist.colname]
    if(physical_dist.is.Mb){ bp <- bp*1000000 }
    save(bp, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/bp.RData'))
    save(bp, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/bp.RData'))
    save(bp, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/bp.RData'))
    
    ## map
    map <- map.this.chr[, map.genetic_dist.colname]
    save(map, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/map.RData'))
    save(map, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/map.RData'))
    save(map, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/map.RData'))
    
    ## subjects
    save(subjects, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/subjects.RData'))
    
    ## strains
    save(strains, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/strains.RData'))
  }
}

#' @export
convert.full.DOQTL.array.to.HAPPY <- function(DOQTL.array, map,
                                              map.locus_name.colname="SNP_ID", map.chr.colname="Chr", map.physical_dist.colname="Mb_NCBI38", map.genetic_dist.colname="cM",
                                              HAPPY.output.path,
                                              bp.column="Mb_NCBI38", physical_dist.is.Mb=TRUE,
                                              allele.labels=LETTERS[1:8],
                                              chr=c(1:19, "X")){
  
  samples <- dimnames(DOQTL.array)[[1]]
  diplotypes <- dimnames(DOQTL.array)[[2]]
  loci <- dimnames(DOQTL.array)[[3]]
  
  full.to.add.matrix <- straineff.mapping.matrix()
  
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/additive/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/full/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  sapply(1:length(chr), function(x) dir.create(path=paste0(HAPPY.output.path, "/genotype/chr", chr[x], "/data/"),
                                               recursive=TRUE, showWarnings=FALSE))
  
  #-------------------------------
  # Marker info
  #-------------------------------
  total.map <- map
  
  for(i in 1:length(loci)){
    chr.locus <- as.character(map[total.map$SNP_ID == loci[i], map.chr.colname])
    
    locus.matrix <- DOQTL.array[,c(1,9,16,22,27,31,34,36,2,3,10,4,11,
                                   17,5,12,18,23,6,13,19,24,28,7,14,
                                   20,25,29,32,8,15,21,26,30,33,35),
                                i]
    
    var_name <- loci[i]
    assign(var_name, locus.matrix)
    full.fn <- paste0(HAPPY.output.path, "/full/chr", chr.locus, "/data/", var_name, ".RData")
    save(list=var_name, file=full.fn)
    
    dosage.matrix <- locus.matrix %*% full.to.add.matrix
    colnames(dosage.matrix) <- allele.labels
    rownames(dosage.matrix) <- rownames(locus.matrix)
    assign(var_name, dosage.matrix)
    add.fn <- paste0(HAPPY.output.path, "/additive/chr", chr.locus, "/data/", var_name, ".RData")
    save(list=var_name, file=add.fn)
  }
  
  subjects <- samples
  strains <- allele.labels
  for(this.chr in chr){
    ## chr
    map.this.chr <- total.map[as.character(total.map[, map.chr.colname]) == this.chr,]
    chromosome <- as.character(map.this.chr$Chr)
    save(chromosome, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/chromosome.RData'))
    save(chromosome, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/chromosome.RData'))
    save(chromosome, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/chromosome.RData'))
    
    ## marker
    markers <- as.character(map.this.chr[, map.locus_name.colname])
    save(markers, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/markers.RData'))
    save(markers, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/markers.RData'))
    save(markers, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/markers.RData'))
    
    ## bp
    bp <- map.this.chr[, map.physical_dist.colname]
    if(physical_dist.is.Mb){ bp <- bp*1000000 }
    save(bp, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/bp.RData'))
    save(bp, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/bp.RData'))
    save(bp, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/bp.RData'))
    
    ## map
    map <- map.this.chr[, map.genetic_dist.colname]
    save(map, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/map.RData'))
    save(map, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/map.RData'))
    save(map, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/map.RData'))
    
    ## subjects
    save(subjects, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/subjects.RData'))
    
    ## strains
    save(strains, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/strains.RData'))
  }
}




