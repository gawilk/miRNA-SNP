# restricts linear model to only TUMOR samples
# 
# IN ADDITION removes samples with: 
#    1) no appreciable gene or miRNA expression
#    2) significant outliers as measured by cook's distance
#
# change to project repo if necessary
load("data/BRCAdata.RData")
load("data/BRCAanova.RData")
source("source/srcFun.R")

#==============================================================================
# filter out unqualified SNPs
#==============================================================================

# remove SNPs with less than 5% genotype representation within TUMORS
# these can deviate from HW and outliers will highly influence regression
classes <- substr(colnames(BRCAdata$gene), 14, 15)
genofreq <- sapply(BRCAdata$mut, function(S) {
  S <- S[classes == "01"]
  tab <- table(S) / sum(table(S))
  all(tab >= 0.05)
})

#get SNPids
snpIDs <- names(genofreq[genofreq])
BRCAanovaTUMOR <- BRCAanova[BRCAanova$snp %in% snpIDs, ]

#==============================================================================
# run on all gene~miRNA:SNP trios
#==============================================================================

# remove samples with no expression
removedExpression <- mcmapply(function(G, M, S) {
  removeOutliers(G = G, M = M, S = S, 
                 data = BRCAdata, method = "expression",
                 threshold = 1, phenotype = "01")
}, BRCAanovaTUMOR$gene, BRCAanovaTUMOR$miRNA, BRCAanovaTUMOR$snp, 
mc.cores = 8, SIMPLIFY = FALSE, USE.NAMES = FALSE)

# remove samples with cookD >= 1
removedCook <- mcmapply(function(G, M, S) {
  removeOutliers(G = G, M = M, S = S, 
                 data = BRCAdata, method = "cook.value", 
                 threshold = 1, phenotype = "01")
}, BRCAanovaTUMOR$gene, BRCAanovaTUMOR$miRNA, BRCAanovaTUMOR$snp,
mc.cores = 8, SIMPLIFY = FALSE, USE.NAMES = FALSE)

#==============================================================================
# p-values and number of outliers
#==============================================================================

# add expression
BRCAanovaTUMOR$pexpr <- sapply(removedExpression, function(X) X$p)
BRCAanovaTUMOR$outliersEXPR <- sapply(removedExpression, function(X) X$outliers)

# add cook
BRCAanovaTUMOR$pcookONE <- sapply(removedCook, function(X) X$p)
BRCAanovaTUMOR$outliersONE <- sapply(removedCook, function(X) X$outliers)

#==============================================================================
# FDR-adjust p-vals
#==============================================================================

BRCAanovaTUMOR$pfdrONE <- p.adjust(BRCAanovaTUMOR$pcookONE, method = "fdr")
BRCAanovaTUMOR$pfdrexpr <- p.adjust(BRCAanovaTUMOR$pexpr, method = "fdr")

#==============================================================================
# save
#==============================================================================

save(BRCAanovaTUMOR, file = "data/BRCAanovaTUMOR.RData")
