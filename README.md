miqtl
====

This package provides a number of different QTL mapping utilities for multiparent populations, both for haplotype association and SNP association.

## Required genetic data format

The functions in this package generally expect that the haplotype data are stored in a genome cache, the directory format output from the HAPPY HMM (Mott *et al.* 2000). The genome cache is a memory efficient approach to storing the genetic data because only a single locus needs to be in memory at a time, rather than a large 3D array, as is output by the HMM of DOQTL and qtl2, and used for subsequent analyses. Currently miqtl cannot use their raw inputs, but does include functions to convert them to a genome cache, and can thus be used for miqtl analyses:

1. convert.DOQTL.to.HAPPY()
2. convert.qtl2.to.HAPPY()

## Reducing genome cache

It is possible that the genome cache is unnecessarily large, with adjacent markers/loci possessing essentially identical information. miqtl includes collapse.genomecache() to create a reduced genome cache, which collapses loci through averaging adjacent markers that are identical according to a set criterion and tolerance level (based on most likely allele or Euclidean distance).
