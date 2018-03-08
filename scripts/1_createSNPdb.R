# get annotated info for affy snp ids & create R db object

# change to project repo if necessary
load("scratch/data/SNPids.RData")
library(pd.genomewidesnp.6)

#==============================================================================
# get info from genomewidesnp
#==============================================================================

get(objects("package:pd.genomewidesnp.6"))
con6 = get(objects("package:pd.genomewidesnp.6"))@getdb()
SNPdb.df <- dbGetQuery(con6, 
	paste0("select * from featureSet where man_fsetid in ('", 
                        paste0(SNPids, collapse = "','"), "')"))

# parse string to get ensembl id and gene symbol: X[1] ensembl, X[5] gene name
geneIDs <- sapply(strsplit(SNPdb.df$gene_assoc, " // "), function(X) {
	c(X[1], X[5])
	})
SNPdb.df$gene <- geneIDs[2, ]
SNPdb.df$ensembl <- geneIDs[1, ]

#==============================================================================
# save database objects
#==============================================================================

# remove extraneous columns and save
SNPdb.df <- subset(SNPdb.df, select = -c(cnv, gene_assoc))
rownames(SNPdb.df) <- SNPdb.df$man_fsetid
save(SNPdb.df, file = "data/SNPdb.df.RData")

# get all SNPs within each gene
SNPbyGene <- split(SNPdb.df[, c("man_fsetid")], SNPdb.df$gene)
save(SNPbyGene, file = "data/SNPbyGene.RData")


