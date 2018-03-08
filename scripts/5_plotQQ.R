# plots qqplot of all SNP:miRNA interactions in TUMORS
#
# change to project repo if necessary
load("data/BRCAanovaTUMOR.RData")
source("source/srcFun.R")

#==============================================================================
# generate plot
#==============================================================================

qq.pretty(na.omit(BRCAanovaTUMOR$pcookONE), "Breast")