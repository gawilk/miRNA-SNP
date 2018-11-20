# runs lm(gene ~ PC + miRNA*SNP) for all gene:miRNA:SNP trios
# PCs are included to account for population substructure
# includes both tumor and normal tissue
#
# change to project repo if necessary
load("data/BRCAanova.RData")
load("data/BRCAdata.RData")
load("data/BRCAevv.RData")
source("source/srcFns.R")
library(parallel)
library(dplyr)

#==============================================================================
# run on trios
#==============================================================================

# run on all combinations
BRCAanovaLMsPCA <- mcmapply(function(G, M, S) {
  runLMPCs(gene = G, PCs = BRCAevv$vectors[,1:2], 
  	mir = M, snp = S, data = BRCAdata)
}, BRCAanova$gene, BRCAanova$miRNA, BRCAanova$snp, 
mc.cores = 16, mc.preschedule = TRUE, USE.NAMES = FALSE, SIMPLIFY = FALSE)

#==============================================================================
# save
#==============================================================================

save(BRCAanovaLMsPCA, file = "data/BRCAanovaLMsPCA.RData")