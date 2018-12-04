# plots regQTL interaction trios w/ PCA substructure correction
# regQTL interactions are plotted for breast tumors
# 
# change to project repo if necessary
load("data/BRCAdata.RData")
load("data/BRCAanovaTUMOR_PCA.RData")
load("data/BRCAevv.RData")
load("data/SNPdb.df.RData")
source("source/srcFns.R")
library(ggplot2)
library(gridExtra)

#==============================================================================
# reorder regQTLs by significance
#==============================================================================

BRCAanovaTUMOR_PCA.ord <- BRCAanovaTUMOR_PCA[order(BRCAanovaTUMOR_PCA$pcookONE), ]

#==============================================================================
# some significant regQTLs
#==============================================================================

g1 <- plotSNPTUMORPCA(gene = BRCAanovaTUMOR_PCA.ord$gene[2], 
                      PCs = BRCAevv$vectors, 
                      mir = BRCAanovaTUMOR_PCA.ord$miRNA[2],
                      snp = BRCAanovaTUMOR_PCA.ord$snp[2], 
                      data = BRCAdata, 
                      method = "cook.value",
                      maintitle = "", 
                      db = SNPdb.df)

g2 <- plotSNPTUMORPCA(gene = BRCAanovaTUMOR_PCA.ord$gene[9], 
                      PCs = BRCAevv$vectors, 
                      mir = BRCAanovaTUMOR_PCA.ord$miRNA[9],
                      snp = BRCAanovaTUMOR_PCA.ord$snp[9], 
                      data = BRCAdata, 
                      method = "cook.value",
                      maintitle = "", 
                      db = SNPdb.df)

g3 <- plotSNPTUMORPCA(gene = BRCAanovaTUMOR_PCA.ord$gene[10], 
                      PCs = BRCAevv$vectors, 
                      mir = BRCAanovaTUMOR_PCA.ord$miRNA[10],
                      snp = BRCAanovaTUMOR_PCA.ord$snp[10], 
                      data = BRCAdata, 
                      method = "cook.value",
                      maintitle = "", 
                      db = SNPdb.df)

BRplots <- list(g1, g2, g3)

#==============================================================================
# plot regQTLs
#==============================================================================

do.call(grid.arrange, c(BRplots, ncol = 3))
