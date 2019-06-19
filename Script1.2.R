library(DESeq2)
library(tidyverse)
library("pasilla")

## Import Data
## -----------

counts <- read.delim("counts.tsv", stringsAsFactors=F)

countdata <- counts %>%
  column_to_rownames("X") %>% # turn the geneid column into rownames
  as.matrix()
head(countdata)

# Get the total number of counts per gene
cpg <- rowSums(countdata, na.rm=F)

# Get the indices of the non expressed genes
nullIndices <- which(cpg==0)

#  Create New matrix with only expressed genes
onlyexpressed <- countdata[-nullIndices,]



## Get meta-data of each sample
## ----------------------------


coldata = read.delim("counts.tsv", row.names = 1)


# separated = strsplit(colnames(onlyexpressed), "_", fixed = TRUE, useBytes = F)
# coldata = do.call(rbind,separated)

# v1 = lapply(separated, '[',1)
# v2 = lapply(separated, '[',2)


# DesEa2 DataObject creation from counts.tsv


dds <- DESeqDataSetFromMatrix (countData = onlyexpressed, # matrix for DESEQ2 analysis
                              colData = coldata, # Rows of colData correspond to columns of countData
                              design = ~ ) #  formula corresponding to each DDSObject, 
# e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment,
#the last variable in the formula for building results tables and plotting

dds