# SNP functions to source for analysis

#==============================================================================
# get all pathway genes from KEGG
#==============================================================================

getPathwayGenes <- function(pathID = NULL) {
  # KEGG pathway genes
  # 
  # Args:
  #  pathID: pathway ID(s), if NULL then get all genes in all pathways
  # Returns:
  #   a list of pathway(s), each entry is vector of genes in pathway
  require(Biobase)
  require(org.Hs.eg.db)
  require(annotate)
  require(KEGG.db)
  keggpath <- org.Hs.egPATH2EG
  keggpath.mapped <- mappedkeys(keggpath)
  paths <- as.list(keggpath[keggpath.mapped])
  if (is.null(pathID)) {
    pathways <- sapply(paths, function(X) {
      eg2symb <- getSYMBOL(X, "org.Hs.eg")
      symbols <- as.character(eg2symb)                       
    })
    return(pathways)
  } else {
    pathways <- sapply(as.character(pathID), function(X) {
      eg2symb <- getSYMBOL(paths[[X]], "org.Hs.eg")
      symbols <- as.character(eg2symb)
    }, simplify = FALSE)
    return(pathways)
  }
}

#==============================================================================
# read SNP data
#==============================================================================

readSNPdata <- function(infile) {
  # reads processed TCGA SNP file into R
  # 
  # Args: 
  #   infile: path of SNP data txt file (after python processing)
  # Returns: 
  #   data frame, SNP IDs by TCGA samples
  SNPdata <- read.table(infile, header = TRUE, row.names = 1, sep = "\t", 
                         stringsAsFactors = FALSE, nrows = 1000000, 
                         check.names = FALSE)
  return(SNPdata)
}

#==============================================================================
# get matching samples for all data
#==============================================================================

getMatchingSamples <- function(mutMat, geneMat, mirMat, barcode.len = 15) {
  # matches TCGA samples, makes list of data objects 
  #
  # Args:
  #   mutMat: matrix of SNP muts, SNP by samples
  #   geneMat: matrix of gene exp, gene by samples
  #   mirMat: matrix of miRNA exp, miRNA by samples
  #   barcode.len: length of barcode to match, default 1-15
  # Returns:
  #   list of matrices with matching samples
  colnames(mutMat) <- substr(colnames(mutMat), 1, barcode.len)
  colnames(geneMat) <- substr(colnames(geneMat), 1, barcode.len)
  colnames(mirMat) <- substr(colnames(mirMat), 1, barcode.len)
  sampnames <- list(muts = colnames(mutMat), 
                    genes = colnames(geneMat), 
                    mirs = colnames(mirMat))
  commonSamps <- Reduce(intersect, sampnames)
  mutMat <- as.matrix(mutMat)[, commonSamps]
  geneMat <- as.matrix(geneMat)[, commonSamps]
  mirMat <- as.matrix(mirMat)[, commonSamps]
  return(list(mut = mutMat, gene = geneMat, mir = mirMat))
}

#==============================================================================
# linear model functions (gene ~ miRNA*SNP)
#==============================================================================

makeDF <- function(miRNA, pathway, pathwayGenes, db) {
  # makes data frame of all gene-mir-snp combinations
  #
  # Args:
  #   miRNA: miRNA name
  #   pathway: pathway KEGG ID
  #   pathwayGenes: pathway list of genes assayed
  #   db: SNPdb split by gene
  # Returns:
  #   df of snp combinations
  Genes <- pathwayGenes[[pathway]]
  SNPs <- db[Genes]
  N <- sapply(SNPs, length)
  data.frame("miRNA" = rep(miRNA, sum(N)), 
             "pathway" = rep(pathway, sum(N)), 
             "gene" = unlist(mapply(rep, Genes, N, SIMPLIFY = FALSE)), 
             "snp" = unlist(SNPs), stringsAsFactors = FALSE)
}

runANOVA <- function(gene, mir, snp, data) {
  # runs linear model, gene ~ miRNA*SNP, for gene-mir-snp combo
  # 
  # Args: 
  #   gene: gene name
  #   mir: mir name
  #   snp: snp id
  #   data: getMatchingSamples output, list of assay data
  # Returns:
  #   interaction p-value of snp modulating mir-gene relationship
  lmfit <- lm(data$gene[gene, ] ~ data$mir[mir, ] * data$mut[[snp]])
  return(anova(lmfit)[3, 5])
}

#==============================================================================
# interaction changes with removed outliers?
#==============================================================================

removeOutliers <- function(G, M, S, data, 
                           method = c("cook.value", "expression"), 
                           threshold = 1, phenotype = "01") {
  # removes samples with no gene or miRNA expression
  # can additionally remove samples using cook's distance
  # recomputes anova interaction p-value
  #
  # Args:
  #  G: gene name (string)
  #  M: mirna name (string)
  #  S: snp ID (affy ID string)
  #  data: data to subset vectors in
  #  method: remove outliers by
  #       "expression" (just samples with no expression) 
  #       "cook.value" (samples with no expression, plus 
  #                   those have large cook's distance)
  #  threshold: threshold for cook's distance, default 1
  #  phenotype: which phenotype to calculate on? default "01" (tumor)
  # Returns:
  #   interaction p-value after outliers are removed
  DF <- data.frame("gene" = data$gene[G, ], 
                   "mir" = data$mir[M, ], 
                   "snp" = data$mut[[S]], stringsAsFactors = FALSE)
  DF <- DF[substr(rownames(DF), 14, 15) == phenotype, ]
  # remove samples with no expression
  minG <- with(data = DF, which(gene <= -13)) #very low gene expression
  minM <- with(data = DF, which(mir == -2)) #no detectable miRNA expression
  noexprsamps <- rownames(DF)[unique(c(minG, minM))]
  if (length(noexprsamps) > 0) {
    DF <- DF[!(rownames(DF) %in% noexprsamps), ]
  }
  # use tryCatch since removing samples may reduce number of genotypes
  # reduction in genotype classes will break anova model 
  if (method == "cook.value") {
    removeSamples <- tryCatch({
      lmfit <- lm(gene ~ mir * snp, data = DF)
      cooksD <- cooks.distance(lmfit)
      outsamps <- names(cooksD[cooksD > threshold])
      alloutliers <- unique(c(noexprsamps, outsamps))
      lmNEW <- lm(gene ~ mir * snp, 
        data = DF[!(rownames(DF) %in% alloutliers), ])
      return(list(p = anova(lmNEW)[3, 5], outliers = length(alloutliers)))
    }, warning = function(war) {
      print(paste("MY WARNING: ", war))
      lmfit <- lm(gene ~ mir * snp, data = DF)
      cooksD <- cooks.distance(lmfit)
      outsamps <- names(cooksD[cooksD > threshold])
      alloutliers <- unique(c(noexprsamps, outsamps))
      lmNEW <- lm(gene ~ mir * snp, 
        data = DF[!(rownames(DF) %in% alloutliers), ])
      return(list(p = anova(lmNEW)[3, 5], outliers = length(alloutliers)))
    }, error = function(err) {
      print(paste("MY ERROR: ", err))
      list("p" = NA, "outliers" = NA)
    }, finally = {})
  } else {
    removeSamples <- tryCatch({
      lmfit <- lm(gene ~ mir * snp, data = DF)
      return(list("p" = anova(lmfit)[3, 5], "outliers" = length(noexprsamps)))
    }, warning = function(war) {
      print(paste("MY WARNING: ", war))
      lmfit <- lm(gene ~ mir * snp, data = DF)
      return(list("p" = anova(lmfit)[3, 5], "outliers" = length(noexprsamps)))
    }, error = function(err) {
      print(paste("MY ERROR: ", err))
      list("p" = NA, "outliers" = NA)
    }, finally = {})
  }
  removeSamples
}

#==============================================================================
# trio plotting functions
#==============================================================================

cleanDF <- function(DF, method = c("expression", "cook.value"),
                    threshold = 1) {
  # remove samples w/o expression or outliers before plotting
  #
  # Args:
  #   DF: data frame 
  #   method: remove outliers by
  #       "expression" (just samples with no expression) 
  #       "cook.value" (samples with no expression, plus 
  #                   those have large cook's distance)
  #   threshold: threshold for cook's distance, default 1
  #  Returns:
  #   cleaned data frame w/ no outlying samples ready for plotting
  minG <- with(data = DF, which(gene <= -13)) #very low gene expression
  minM <- with(data = DF, which(mir == -2)) #no detectable miRNA expression
  noexprsamps <- rownames(DF)[unique(c(minG, minM))]
  if (length(noexprsamps) > 0) {
    DF <- DF[!(rownames(DF) %in% noexprsamps), ]
  }
  if (method == "expression") {
    return(DF)
  } else {
    lmfit <- lm(gene ~ mir * mut, data = DF)
    cooksD <- cooks.distance(lmfit)
    outsamps <- names(cooksD[cooksD > threshold])
    DF <- DF[!(rownames(DF) %in% outsamps), ]
    return(DF)
  }
}

plotSNPTUMOR <- function(miRNA, gene, mut, data, filter = TRUE, 
                         method = c("expression", "cook.value"), 
                         threshold = 1, maintitle, db) {
  # plots SNP interaction plot for only tumor samples
  # changes SNP Affy ID to rsID
  # also changes 0,1,2 encoding to actual genotypes
  #
  # Args:
  #  miRNA: miRNA ID string 
  #  gene: gene name string
  #  mut: SNP Affy ID
  #  data: output from getMatchingSamples
  #  filter: remove outliers? 
  #  method: choose which method
  #  threshold: value of cook's threshold
  #  maintitle: title of plot
  #  db: SNPdb.df for converting SNP Affy ID to rsID
  # Returns:
  #  trio interaction plot 
  require(ggplot2)
  require(gridExtra)
  #prep dataframe for plotting
  DF <- data.frame(mir = data$mir[miRNA, ], 
                   gene = data$gene[gene, ],
                   mut = data$mut[[mut]], 
                   stringsAsFactors = FALSE)
  DF <- DF[substr(rownames(DF), 14, 15) == "01", ]
  DF <- DF[!is.na(DF$mut), ]
  DF <- cleanDF(DF, method = method, threshold = threshold)
  #get SNP genotype and rsID
  snpinfo <- db[db$man_fsetid == mut, ]
  genotypes <- c(paste0(rep(snpinfo$allele_a, 2), collapse = ""), 
                 paste0(snpinfo$allele_a, snpinfo$allele_b), 
                 paste0(rep(snpinfo$allele_b, 2), collapse = ""))
  names(genotypes) <- c("0", "1", "2")
  DF$mut <- as.factor(genotypes[as.character(DF$mut)])
  genoTable <- table(DF$mut)
  g <- ggplot(DF, aes(x = mir, y = gene)) +
    geom_point(aes(shape = mut, color = mut), size = 2, alpha = 1) + 
    scale_shape(solid = FALSE) + 
    geom_line(aes(color = mut, linetype = mut), stat = "smooth", 
              alpha = 0.5, method = "lm", 
              show.legend = FALSE, se = FALSE, size = 1) + 
    scale_color_manual(as.character(snpinfo$dbsnp_rs_id), 
                       values = c("red", "blue", "forestgreen"), 
                       labels = mapply(function(X, Y) {
                         paste0(X, " (", Y, ")")
                       }, names(genoTable), as.vector(genoTable))) + 
    scale_shape_manual(as.character(snpinfo$dbsnp_rs_id), 
                       values = c(0, 1, 2), 
                       labels = mapply(function(X, Y) {
                         paste0(X, " (", Y, ")")
                       }, names(genoTable), as.vector(genoTable))) + 
    scale_linetype_manual(as.character(snpinfo$dbsnp_rs_id), 
                          values = c(1, 1, 1), 
                          labels = mapply(function(X, Y) {
                            paste0(X, " (", Y, ")")
                          }, names(genoTable), as.vector(genoTable))) +
    xlab(miRNA) + ylab(gene) +  
    ggtitle(paste0(as.character(maintitle))) + 
    theme(panel.background = element_rect(fill = "white", colour = "white"), 
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16), 
          title = element_text(size = 14),
          axis.line.x = element_line(size = 0.5),
          axis.line.y = element_line(size = 0.5), 
          plot.title = element_text(hjust = 0.5))
  return(g)
}

#==============================================================================
# SNP in Hardy-Weinberg Equilibrium?
#==============================================================================

SNPinHW <- function(snp) {
  # test if SNP is in HW
  # 
  # Args:
  #  snp: snp vector
  # Returns:
  #  is particular snp in hardy-weinberg equilibrium?
  #  p-value of goodness of fit test
  #  ordered from "0" to "2" SNP genotype encoding
  snp <- as.integer(as.character(na.omit(snp)))
  O <- as.vector(table(snp)) #ordered from 0 to 2
  p <- mean(snp) / 2 #major allele frequency
  q <- 1 - p #minor allele frequency
  E <- c(q^2, 2*p*q, p^2) * length(snp)
  if (length(O) != length(E)) {
    warning("not all genotypes represented in snp")
    return(list("chisq" = NA, "p" = NA))
  }
  stat <- sum((O - E)^2 / E)
  pval <- pchisq(q = stat, df = 1, lower.tail = FALSE)
  return(list("chisq" = stat, "p" = pval))
}

#==============================================================================
# chromosome absolute positions
#==============================================================================

getCHRcoords <- function(df) {
  # get absolute chromosome coords
  #
  # Args: 
  #   df: full data frame with SNP range
  # Returns: 
  #   absolute coords of chromosomes in genome
  chrs <- paste0("chr", c(as.character(seq(1, 22)), "X", "Y"))
  DFsplit <- split(df, as.factor(df$chrom))
  sizes <- sapply(chrs, function(X) max(DFsplit[[X]]$pos))
  sizes <- sizes[is.finite(sizes)]
  sums <- cumsum(as.numeric(sizes))
  starts <- c(0, sums[(1:length(sums)) - 1])
  chrmap <- data.frame(chr = paste0("chr", c(as.character(seq(1, 22)), "X")), 
                       start = starts)
  rownames(chrmap) <- chrmap$chr
  return(chrmap)
}

#==============================================================================
# chromosome absolute positions
#==============================================================================

qq.pretty <- function(pvals, maintitle) {
  # plot distribution of interaction p-values
  #
  # Args:
  #   pvals: vector of pvals (NAs omitted)
  #   maintitle: title of plot
  # Returns:
  #   qqplot w/ genome-wide bonferroni correction line
  obs.p <- sort(pvals)
  exp.p <- seq_along(pvals) / (length(pvals) + 1)
  plot(-log10(exp.p), -log10(obs.p), 
       main = maintitle, xlab = "", ylab = "", 
       cex = 2, cex.main = 2, cex.axis = 1.75)
  abline(a = 0, b = 1, col = 2, lwd = 2)
  abline(a = -log10(0.05/length(pvals)), b = 0, 
         col = "blue", lwd = 2)
  mtext(expression("Expected "*-log[10]*"p"), side = 1, 
    line = 3, cex = 1.75)
  mtext(expression("Observed "*-log[10]*"p"), side = 2, 
    line = 2.5, cex = 1.75)
}

#==============================================================================
# convert SNP data into matrix form
#==============================================================================

convertToMatrix <- function(InputData) {
  # converts SNP data into matrix
  #
  # Args:
  #   InputData: list of objects with common tumor samples
  # Returns:
  #   
  DF <- as.data.frame(lapply(InputData$mut, function(X) {
    as.numeric(as.character(X))
  }), stringsAsFactors = FALSE, optional = TRUE)
  rownames(DF) <- colnames(InputData$gene)
  mat <- as.matrix(DF)
  return(mat)
}