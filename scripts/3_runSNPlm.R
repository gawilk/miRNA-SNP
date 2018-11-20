# finds unique gene:mir:SNP combinations on dysregulated KEGG pathways
# this script references breast cancer data
#      
# script assumes miRNA-pathway dysregulation stats are already available 
#     1) "BRCAstats1e5.RData" dataset contains pathway dysregulation by 
#        miRNAs in tumors
# data object can be generated using our miRNA-pathways repo:
# github: https://github.com/gawilk/miRNA-pathways
#
# change to project repo if necessary
load("data/SNPbyGene.RData")
load("data/BRCAstats1e5.RData")
source("source/srcFns.R")
library(parallel)
library(dplyr)

#==============================================================================
# clean objects & get pairs
#==============================================================================

# filters genes & SNPs to only those assayed in data
pathGenes <- lapply(getPathwayGenes(), function(P) {
	intersect(P, rownames(BRCAdata$gene))
})
SNPbyGene.uni <- SNPbyGene[unique(unlist(pathGenes))]
SNPbyGene.uni.int <- lapply(SNPbyGene.uni, function(S) {
	intersect(S, names(BRCAdata$mut))
})
SNPbyGene.filt <- SNPbyGene.uni.int[sapply(SNPbyGene.uni.int, length) > 0]
pathGenes.filt <- lapply(pathGenes, function(P) {
	intersect(P, names(SNPbyGene.filt))
})

# get most significant miR-path pairs
BRCAstats1e5.sig <- subset(BRCAstats1e5, pval < 0.01)

#==============================================================================
# precompute all unique gene:mir:snp combinations
#==============================================================================

# df of all combinations to run
BRCAanova <- mapply(function(M, P) {
  makeDF(miRNA = M, pathway = P, 
  	pathwayGenes = pathGenes.filt, db = SNPbyGene.filt)
}, BRCAstats1e5.sig$miRNA, BRCAstats1e5.sig$pathway, SIMPLIFY = FALSE)
BRCAanova <- rbind_all(BRCAanova)
BRCAanova <- BRCAanova[!duplicated(BRCAanova[, c("miRNA", "gene", "snp")]), ]

# garbage collect
rm(BRCAstats1e5, SNPbyGene, SNPbyGene.uni, SNPbyGene.uni.int, pathGenes)

#==============================================================================
# save combos
#==============================================================================

# save
save(BRCAanova, file = "data/BRCAanova.RData")

