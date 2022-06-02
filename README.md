miQTL
====

This package provides a number of different QTL mapping utilities for multiparent populations, both for haplotype association and SNP association.

## Required genetic data format

The functions in this package generally expect that the haplotype data are stored in a genome cache, the directory format output from the HAPPY HMM (Mott *et al.* 2000). The genome cache is a memory efficient approach to storing the genetic data because only a single locus needs to be in memory at a time, rather than a large 3D array, as is output by the HMM of DOQTL and qtl2, and used for subsequent analyses. Currently miqtl cannot use their raw inputs, but does include functions to convert them to a genome cache, and can thus be used for miqtl analyses:

1. convert.DOQTL.to.HAPPY()
2. convert.qtl2.to.HAPPY()

## Reducing genome cache

It is possible that the genome cache is unnecessarily large, with adjacent markers/loci possessing essentially identical information. miqtl includes collapse.genomecache() to create a reduced genome cache, which collapses loci through averaging adjacent markers that are identical according to a set criterion and tolerance level (based on most likely allele or Euclidean distance).

## Installation

```r

#install.packages("remotes")

remotes::install_github("gkeele/miqtl")

```

## Simulating Data

We use the package [sparc](https://github.com/gkeele/sparcc) to simulate phenotypes from a multiparent population, [Collaborative Cross](https://pubmed.ncbi.nlm.nih.gov/21411855/) (CC).This also requires the [sparc cache repository](https://github.com/gkeele/sparcc_cache) to ensure the simulated data is representative of the Collaborative Cross strains. 

```r

#devtools::install_github("gkeele/sparcc")

library(sparcc)

## Download the sparc genome cache from GitHub
download.file(url = "https://github.com/gkeele/sparcc_cache/archive/master.zip", destfile = "~/sparcc_cache.zip")
unzip(zipfile = "~/sparcc_cache.zip")
genomecache <- "~/sparcc_cache"

phenotypes <- sim.CC.data(genomecache = genomecache, 
                          num.lines = 40, 
                          num.sim = 10, 
                          num.replicates = 2, 
                          qtl.effect.size = 0.8)

```


## Example workflow

The following code demonstrates a simple example workflow to implement miQTL. This gives a demonstration of a QTL scan in the Collaborative Cross population and compares scans which handle haplotype uncertainty using regression on probabilities (ROP) and multiple imputation. 

```r
# QTL scan using ROP

miqtl.rop.scan <- scan.h2lmm(genomecache = genomecache, 
                             data = phenotypes$data, 
                             pheno.id = "SUBJECT.NAME.1", 
                             geno.id =  "SUBJECT.NAME.1", 
                             formula = sim.y.1 ~ 1, 
                             use.multi.impute = FALSE, 
                             return.allele.effects = TRUE)

genome.plotter.whole(scan.list=list(ROP = miqtl.rop.scan))

# QTL scan using multiple imputation

miqtl.mi11.scan <- scan.h2lmm(genomecache = genomecache, 
                              data = phenotypes$data, 
                              pheno.id = "SUBJECT.NAME.1", 
                              geno.id = "SUBJECT.NAME.1", 
                              formula=sim.y.1 ~ 1, 
                              use.multi.impute = TRUE, 
                              num.imp = 11)

genome.plotter.whole(scan.list = list(miqtl.mi11.scan),i)

```
miQTL generates significance thresholds for each scan using a permutation test

```r
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = miqtl.rop.scan, 
                                                      method = "permutation", num.samples = 10)
                                                      
permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes$data,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1)

permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, 
                                        percentile = 0.95)
                                        
genome.plotter.whole(scan.list=list(ROP = miqtl.rop.scan), 
                     hard.thresholds = permute_threshold)

```
The top peak can be isolated using grab.locus.from.scan() be returning the locus with the minimum peak. The positional confidence of that QTL can be identified using parametric bootstrap analysis. 

```r

qtl_locus <- grab.locus.from.scan(miqtl.rop.scan)

## Heatmap of QTL 

prob.heatmap(marker = qtl_locus, 
             p.value = min(miqtl.rop.scan$p.value), 
             genomecache = genomecache, 
             model = "additive", 
             phenotype = "sim.y.1", 
             phenotype.data = phenotypes$data %>% dplyr::rename(SUBJECT.NAME = SUBJECT.NAME.1))
             
## QTL positional confidence interval using parametric bootstraps

qtl.fit <- scan.h2lmm(genomecache = genomecache, 
                      data = phenotypes$data, 
                      pheno.id = "SUBJECT.NAME.1", 
                      geno.id =  "SUBJECT.NAME.1",
                      formula = sim.y.1 ~ 1, 
                      use.multi.impute=FALSE, 
                      just.these.loci = qtl_locus)
                      
qtl.locus.bootstraps <- generate.sample.outcomes.matrix(scan.object = qtl.fit, 
                                                        model.type = "alt", 
                                                        method = "bootstrap", 
                                                        num.sample = 100, 
                                                        seed = 1, 
                                                        use.REML = TRUE, 
                                                        use.BLUP = TRUE)
                                                        
qtl.bootstrap.scans <- run.positional.scans(sim.object = qtl.locus.bootstraps, 
                                            keep.full.scans=TRUE,
                                            genomecache = genomecache, 
                                            data = phenotypes$data,
                                            use.multi.impute = FALSE, 
                                            scan.seed = 1)

single.chr.plotter.w.ci(scan.object = miqtl.rop.scan, 
                        qtl.ci.object = qtl.bootstrap.scans, 
                        ci.type = "bootstrap", 
                        scan.type = "ROP", 
                        scale = "Mb", 
                        these.col = rep("#7BAFD4", 100))

```

Finally, haplotype effects of the region of interest can be visualized to betteer understand founder strain contributions to alleles driving variation.

```r

## QTL haplotype effects in region of QTL

genome.plotter.region(haplotype.association = list(miqtl.rop.scan), 
                      chr = 14)
                      
allele.plotter.region(scan.object = miqtl.rop.scan, chr = 14)


```

