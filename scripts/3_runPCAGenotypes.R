# computes PCA of SNP genotypes to correct for population substructure

# change to project repo if necessary
load("data/BRCAdata.RData")
load("data/BRCAesetTPM.RData")
source("source/srcFns.R")
library(snpStats)

#==============================================================================
# compute PCs
#==============================================================================

mat <- convertToMatrix(BRCAdata)
newMat <- new("SnpMatrix", mat)
xxtmat <- xxt(newMat, correct.for.missing = FALSE)
BRCAevv <- eigen(xxtmat)

#==============================================================================
# save
#==============================================================================

# save decomposition
save(BRCAevv, file = "data/BRCAevv.RData")

# get races of samples & save
BRCAcols <- getRaces(mat, BRCAesetTPM)
save(BRCAcols, file = "data/BRCAcols.RData")
