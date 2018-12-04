# generates QQ-plot of all SNP:miRNA interactions (regQTLs) in TUMOR samples

# change to project repo if necessary
load("data/BRCAanovaTUMOR_PCA.RData")
source("source/srcFun.R")

#==============================================================================
# generate plot
#==============================================================================

qq.pretty(na.omit(BRCAanovaTUMOR_PCA$pcookONE), "Breast")