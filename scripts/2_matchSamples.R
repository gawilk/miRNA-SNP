# finds data for common tumor samples across RNAseq, miRNAseq, and SNP genotypes 
# creates list of data objects to link them together
#
# script assumes RNAseq & miRNAseq for samples already available 
#     1) BRCAesetTPM.RData is the RNAseq dataset 
#     2) BRCA_miRNASeq_HiSeq.Rda is the miRNAseq dataset
# both objects can be generated using our miRNA-pathways repo:
# github: https://github.com/gawilk/miRNA-pathways
#
# change to project repo if necessary
load("data/BRCA.SNP.RData")
load("data/BRCAesetTPM.RData")
load("data/BRCA_miRNASeq_HiSeq.Rda")
load("data/SNPbyGene.RData")
source("source/srcFns.R")
library(Biobase)

#==============================================================================
# get matching samples for all three data types
#==============================================================================

# filter genes
# get only genes assayed with SNP data which are on KEGG pathways
Genes <- Reduce(intersect, list(featureNames(BRCAesetTPM), 
                                names(SNPbyGene), 
                                unique(unlist(getPathwayGenes()))))

# filter muts
# get muts on those filtered genes which have been assayed
Muts <- intersect(rownames(BRCA.SNP), unique(unlist(SNPbyGene[Genes])))

# filter miRNA data 
# remove miRNAs with very low expression & log2 transform
BRCA_miRNA <- BRCA_miRNASeq_HiSeq[apply(BRCA_miRNASeq_HiSeq, 1, function(X) {
	sum(X > 1) > (0.5 * length(X))
}), ]

# combine all
BRCAdata <- getMatchingSamples(mutMat = BRCA.SNP[Muts, ], 
                           geneMat = exprs(BRCAesetTPM[Genes, ]), 
                           mirMat = log2(as.matrix(BRCA_miRNA) + 0.25),
                           barcode.len = 15)

#==============================================================================
# filter SNPs
#==============================================================================

# filter SNPs by minor allele frequency (MAF)
# minor allele should be be 1% or greater
# need to add bounds to both sides!
BRCAdata$mut <- BRCAdata$mut[apply(BRCAdata$mut, 1, function(X) {
  MAF <- mean(X, na.rm = TRUE)/2
  (MAF <= 0.99) & (MAF >= 0.01)
}), ]
mutmaf <- nrow(BRCAdata$mut)

# remove SNPs with less than 3 genotypes
BRCAdata$mut <- BRCAdata$mut[apply(BRCAdata$mut, 1, function(X) {
	length(levels(as.factor(X))) == 3
}), ]

# encode SNPs as factors
BRCAdata$mut <- lapply(as.data.frame(t(BRCAdata$mut)), as.factor)

#==============================================================================
# save object
#==============================================================================

save(BRCAdata, file = "data/BRCAdata.RData")

#==============================================================================
# print filtration info
#==============================================================================

paste("starting with", nrow(BRCA.SNP), "unique SNPs")
paste("SNPs on assayed genes in KEGG pathways reduces to", 
	length(Muts), "SNPs")
paste("1% minor allele frequency filtering reduces to", mutmaf, "SNPs")
paste("removing SNPs with less than 3 genotypes reduces to", 
	length(BRCAdata$mut), "SNPs")
