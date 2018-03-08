# reads BRCA SNP data into R data frame
# this script assumes: 
#     1) SNP genotype data has already been processed by python script
#        and is R-readable
#     2) the processed txt file resides in data/ in the repo
#
# change to project repo if necessary
source("source/srcFns.R")

#==============================================================================
# create & save data objects
#==============================================================================

# read into dataframe
BRCA.SNP <- readSNPdata("data/BRCA_SNP.txt")
save(BRCA.SNP, file = "data/BRCA.SNP.RData")

# save SNP Affy IDs as well
SNPids <- rownames(BRCA.SNP)
save(SNPids, file = "data/SNPids.RData")