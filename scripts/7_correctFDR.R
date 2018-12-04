# adds tumor lm pvals to dataframe with outlier removal (using cook's distance)
# FDR-correct pvals for the large number of tests

# change to project repo if necessary
load("data/BRCAanovaTUMOR.RData")
load("data/BRCAremovedCookPCA.RData")

#==============================================================================
# add cols to dataframe
#==============================================================================

# make copy of dataframe to ensure no overwriting & garbage collect
BRCAanovaTUMOR_PCA <- BRCAanovaTUMOR
rm(BRCAanovaTUMOR)

# lm with tumor 
BRCAanovaTUMOR_PCA$pcookONE <- sapply(BRCAremovedCookPCA, function(X) X$p)
BRCAanovaTUMOR_PCA$outliersONE <- sapply(BRCAremovedCookPCA, function(X) {
	X$outliers
})

# FDR-adjust pvals
BRCAanovaTUMOR_PCA$pfdrONE <- p.adjust(BRCAanovaTUMOR_PCA$pcookONE, 
	method = "fdr")

#==============================================================================
# save
#==============================================================================

save(BRCAanovaTUMOR_PCA, file = "data/BRCAanovaTUMOR_PCA.RData")
