# generates manhattan plot of all gene~miRNA:SNP interactions
#
# change to project repo if necessary
load("data/BRCAanovaTUMOR.RData")
load("data/SNPdb.df.RData")
source("source/srcFun.R")
library(reshape2)
library(ggplot2)

#==============================================================================
# add info to anova df
#==============================================================================

BRCAanovaTUMOR$rsID <- SNPdb.df[BRCAanovaTUMOR$snp, ]$dbsnp_rs_id
BRCAanovaTUMOR$pos <- SNPdb.df[BRCAanovaTUMOR$snp, ]$physical_pos
BRCAanovaTUMOR$chrom <- paste0("chr", SNPdb.df[BRCAanovaTUMOR$snp, ]$chrom)

#==============================================================================
# map absolute positions of SNPs in genome
#==============================================================================

chrcoords <- getCHRcoords(BRCAanovaTUMOR)
BRCAanovaTUMOR$abspos <- BRCAanovaTUMOR$pos + 
		chrcoords[BRCAanovaTUMOR$chrom, ]$start

#==============================================================================
# save changes
#==============================================================================

save(BRCAanovaTUMOR, file = "data/BRCAanovaTUMOR.RData")

#==============================================================================
# generate manhattan plot
#==============================================================================

ggBreast <- ggplot(BRCAanovaTUMOR)
ggBreast <- ggBreast + geom_point(aes(x = abspos, 
	y = -log10(pfdrONE), color = chrom), size = 1.5, alpha = 0.5)
ggBreast <- ggBreast + scale_x_continuous(labels = rownames(chrcoords), 
	breaks = chrcoords$start)
ggBreast <- ggBreast + theme(panel.background = element_rect(fill = "white", 
	colour = "white"), legend.position = "none", 
               axis.text.x = element_text(size = 16, angle = 90, hjust = 1), 
               axis.text.y = element_text(size = 16), 
               axis.title = element_text(size = 16))
ggBreast <- ggBreast + xlab("") + ylab(expression(-log[10](p[FDR])))
ggBreast