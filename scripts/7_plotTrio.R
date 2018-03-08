# plots a gene~miRNA:SNP interaction trio
# 
# change to project repo if necessary
load("data/BRCAanovaTUMOR.RData")
load("data/SNPdb.df.RData")
source("source/srcFun.R")
library(ggplot2)

#==============================================================================
# plot trio
#==============================================================================

# pick any trio! here's trio from 5th row
combo <- BRCAanovaTUMOR[5, ]
g1 <- plotSNPTUMOR(miRNA = combo$miRNA, gene = combo$gene, 
             mut = combo$snp, data = BRCAdata, filter = TRUE, 
             method = "cook.value", threshold = 1, 
             maintitle = "", db = SNPdb.df)
g1