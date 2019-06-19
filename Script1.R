library(DESeq2)
library(tidyverse)


## Import Data
## -----------

counts <- read.delim("counts.tsv", stringsAsFactors=F)

countdata <- counts %>%
  column_to_rownames("X") %>% # turn the geneid column into rownames
  as.matrix()
head(countdata)

# Get the total number of counts per gene
cpg <- rowSums(countdata, na.rm=F)

# Test the nullity of expression of each gene
nullGenes_test <- cpg==0

# Get the indices of the non expressed genes
nullIndices <- which(nullGenes_test)

# Get the total number of non expressed genes
nullNumber = length(nullIndices)

# New matrix with only expressed genes
onlyexpressed <- countdata[-nullIndices,]

# Choose only well expressed genes
# keep <- rowSums(, na.rm=F) >500000
# countdata <- countdata[keep,]
# dim(countdata)


## Get meta-data of each sample
## ----------------------------
separated = strsplit(colnames(onlyexpressed), "_", fixed = TRUE, useBytes = F)
do.call(rbind,separated)


# v1 = lapply(separated, '[',1)
# v2 = lapply(separated, '[',2)
# for (i in v1) separated2 = strsplit(v1[i], " ")



# DesEa2 DataObject creation from counts.tsv

dds <- DESeqDataSetFromMatrix(countData = onlyexpressed, # matrix for DESEQ2 analysis
                              colData = coldata, # Rows of colData correspond to columns of countData
                              design = ~ condition) #  formula corresponding to each DDSObject, 
# e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment,
#the last variable in the formula for building results tables and plotting

dds










#Bar plot of genes expression 
expression <- rowSums(countdata)
barplot (expression,
         names = names(expression),
         las=2,
         main="Barplot of genes expression")
abline(h = 20e6, lty= 2)

# Get log2 counts per million
logcounts <- log2(countdata + 1)
