# generates manhattan plot of all regQTL interactions in breast tumors

# change to project repo if necessary
load("data/BRCAanovaTUMOR_PCA.RData")
load("data/SNPdb.df.RData")
source("source/srcFun.R")
library(reshape2)
library(ggplot2)

#==============================================================================
# add rsIDs, SNP positions, and chromosome data
#==============================================================================

BRCAanovaTUMOR_PCA$rsID <- SNPdb.df[BRCAanovaTUMOR_PCA$snp, ]$dbsnp_rs_id
BRCAanovaTUMOR_PCA$pos <- SNPdb.df[BRCAanovaTUMOR_PCA$snp, ]$physical_pos
BRCAanovaTUMOR_PCA$chrom <- paste0("chr", 
	SNPdb.df[BRCAanovaTUMOR_PCA$snp, ]$chrom)

#==============================================================================
# map SNP positions to absolute positions in genome
#==============================================================================

chrcoords <- getCHRcoords(BRCAanovaTUMOR_PCA)
BRCAanovaTUMOR_PCA$abspos <- BRCAanovaTUMOR_PCA$pos + 
		chrcoords[BRCAanovaTUMOR_PCA$chrom, ]$start

#==============================================================================
# save changes
#==============================================================================

save(BRCAanovaTUMOR_PCA, file = "data/BRCAanovaTUMOR_PCA.RData")

#==============================================================================
# generate manhattan plot
#==============================================================================

# create threshold to magnify significant regQTLs
larger <- with(BRCAanovaTUMOR_PCA, -log10(pfdrONE) > 2)

# plot!
ggBreast <- ggplot(BRCAanovaTUMOR_PCA)
ggBreast <- ggBreast + geom_point(aes(x = abspos, 
	y = -log10(pfdrONE), color = chrom), size = 0.7, alpha = 0.5)
ggBreast <- ggBreast + geom_point(data=BRCAanovaTUMOR_PCA[larger, ], 
                                  aes(x = abspos, y = -log10(pfdrONE), 
                                  	color = chrom), size = 1.5, alpha = 0.8)
ggBreast <- ggBreast + scale_x_continuous(labels = rownames(chrcoords), 
                                          breaks = chrcoords$start)
ggBreast <- ggBreast + 
	theme(panel.background = element_rect(fill = "white", colour = "white"), 
          legend.position = "none",
          axis.line = element_line(colour = "black", size = 1, 
                                   linetype = "solid", lineend = "square"),
          panel.grid.major.y = element_line(colour = "grey50", size = 0.15),
          axis.text.x = element_text(size = 8, angle = 65, hjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title = element_text(size = 8))
ggBreast <- ggBreast + xlab("") + ylab(expression(-log[10](p[FDR])))
ggBreast
