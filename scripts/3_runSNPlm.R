# runs linear model for all unique gene~miRNA:SNP unique combinations 
# found on dysregulated pathways by miRNAs in tumors
#
# NOTE! because of large size, uses multiple cores to run code, can modify
#       to run on a single core
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

# filters genes/SNPs to only those assayed in data
pathGenes <- lapply(getPathwayGenes(), function(P) {
	intersect(P, rownames(BRCAdata$gene))
})
SNPbyGene <- SNPbyGene[unique(unlist(pathGenes))]
SNPbyGene <- lapply(SNPbyGene, function(S) intersect(S, names(BRCAdata$mut)))
SNPbyGene <- SNPbyGene[sapply(SNPbyGene, length) > 0]
pathGenes <- lapply(pathGenes, function(P) intersect(P, names(SNPbyGene)))

# get most significant miR-path pairs
BRCAstats1e5 <- subset(BRCAstats1e5, pval < 0.01)

#==============================================================================
# precompute all unique gene-mir-snp combinations
#==============================================================================

# df of all combinations to run
BRCAanova <- mapply(function(M, P) {
  makeDF(miRNA = M, pathway = P, pathwayGenes = pathGenes, db = SNPbyGene)
}, BRCAstats1e5$miRNA, BRCAstats1e5$pathway, SIMPLIFY = FALSE)
BRCAanova <- rbind_all(BRCAanova)
BRCAanova <- BRCAanova[!duplicated(BRCAanova[, c("miRNA", "gene", "snp")]), ]

# remove objects
rm(BRCAstats1e5, SNPbyGene, pathGenes)

#==============================================================================
# run on significant pairs
#==============================================================================

# run on all combinations
BRCAanova$pval <- mcmapply(function(G, M, S) {
  runANOVA(gene = G, mir = M, snp = S, data = BRCAdata)
}, BRCAanova$gene, BRCAanova$miRNA, 
BRCAanova$snp, mc.cores = 8, USE.NAMES = FALSE)

# save
save(BRCAanova, file = "data/BRCAanova.RData")

