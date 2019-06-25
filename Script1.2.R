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


coldata = colnames(onlyexpressed)

separated = strsplit((coldata), "_", fixed = TRUE, useBytes = F)
coldata = cbind(coldata, do.call(rbind,separated))

# v1 = lapply(separated, '[',1)
# v2 = lapply(separated, '[',2)



## Deseq2 data object created
## ----------------------------
colnames(coldata) [c(2,3)] <- c("id", "day")
dds <- DESeqDataSetFromMatrix (countData = onlyexpressed, # matrix for DESEQ2 analysis
                              colData = coldata, # Rows of colData correspond to columns of countData
                              design = V2 ~ V3 ) #  V2 is subject number anf V2 is day from vaccination

#formula corresponding to each DDSObject, 
# e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment,
#the last variable in the formula for building results tables and plotting


# explicit reference level 
dds$V2 = relevel(dds$V3, ref = "d0" )

# results d0 vs other days :
dds <- DESeq(dds)
res = results(dds)
resultsNames(dds) # different days (d1,d3,d7) comparison vs d0

# result for adjusted p, alpha =0.05 cause deafault is 0.1 alpha
res05 <- results(dds, alpha=0.05)

# remove noise = log fold change shrinkage (always for biological samples )
resLFC = lfcShrink(dds,  coef = "V3_d1_vs_d0", type = "apeglm")
# apeglm method for effect size shrinkage improving on the previous estimator

# order results by smallest p-value 
resOrdered = res[order(res$pvalue),]

# How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=T)



## Independent Hypothesis Weighting
## ----------------------------

#Data-driven hypothesis weighting increases detection power in genome-scale mult testing
library ("IHW")
resIHW = results (dds, filterFun = ihw)
summary (resIHW)


## Exploring and exporting results
## ----------------------------

# MA-plot for alpha 0.05 non shrunk results
plotMA(res05, ylim=c(-2,2))
# MA-plot for shrunken log2 LFC, that removes noise
plotMA(resLFC, ylim=c(-2,2))
# use the function identify to interactively detect the row number of individual genes 
# by clicking on the plot
# idx = identify(res$baseMean, res$log2FoldChange)
# rownames(res) [idx]


# Alternative shrinkage estimators


#------------------ Code Hadrien
counts <- DESeq::counts(dds)
toKeep <- which(SummarizedExperiment::colData(se)$day %in% c(0,1,3,7))
coldataSE <- SummarizedExperiment::colData(se)[toKeep,]
day <- factor(coldataSE$day)#,labels = c("_d0","_d1","_d3","_d7"))
id_num <- factor(coldataSE$id_num)
group <- rep(1,length(id_num))
group[which(as.numeric(levels(id_num)[id_num])>10)] <-2
design <- model.matri~0+day);colnames(design);rownames(design) <-
  colnames(counts)









