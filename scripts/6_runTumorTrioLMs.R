# restricts linear model lm(gene ~ PC + miRNA*SNP) to only TUMOR samples
# 
# IN ADDITION removes samples with: 
#    1) no appreciable gene or miRNA expression
#    2) significant outliers as measured by Cook's distance
#
# change to project repo if necessary
load("data/BRCAdata.RData")
load("data/BRCAanova.RData")
load("data/BRCAevv.RData")
source("source/srcFun.R")
library(parallel)
library(dplyr)

#==============================================================================
# filter out SNPs with low genotype frequencies in tumors
#==============================================================================

# remove SNPs < 5% genotype frequency within TUMORS
# these can deviate from HW and outliers will highly influence regression
classes <- substr(colnames(BRCAdata$gene), 14, 15)
genofreq <- sapply(BRCAdata$mut, function(S) {
  S <- S[classes == "01"]
  tab <- table(S) / sum(table(S))
  all(tab >= 0.05)
})

# retrieve relevant SNPids & subset
snpIDs <- names(genofreq[genofreq])
BRCAanovaTUMOR <- BRCAanova[BRCAanova$snp %in% snpIDs, ]

#==============================================================================
# save filtered trios
#==============================================================================

# save tumor trios & garbage collect full set
save(BRCAanovaTUMOR, file = "data/BRCAanovaTUMOR.RData")
rm(BRCAanova)

#==============================================================================
# run new lm with Cook's distance correction
#==============================================================================

# remove samples with cookD >= 1
BRCAremovedCookPCA <- mcmapply(function(G, M, S) {
  removeOutliersPCA(G = G, M = M, S = S, PCs = BRCAevv$vectors[, 1:2],
                 data = BRCAdata, method = "cook.value", 
                 threshold = 1, phenotype = "01")
}, BRCAanovaTUMOR$gene, BRCAanovaTUMOR$miRNA, BRCAanovaTUMOR$snp,
mc.cores = 16, mc.preschedule = TRUE, SIMPLIFY = FALSE, USE.NAMES = FALSE)

#==============================================================================
# save
#==============================================================================

save(BRCAremovedCookPCA, file = "data/BRCAremovedCookPCA.RData")
