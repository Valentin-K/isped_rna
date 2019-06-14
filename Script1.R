library(DESeq2)
library(tidyverse)
counts <- read.delim("counts.tsv", stringsAsFactors=F)

countdata <- counts %>%
  column_to_rownames("X") %>% # turn the geneid column into rownames
  as.matrix()
head(countdata)


keep <- rowSums(countdata, na.rm=F) >500000
countdata <- countdata[keep,]
dim(countdata)

#Bar plot of genes expression 
expression <- rowSums(countdata)
barplot (expression,
         names = names(expression),
         las=2,
         main="Barplot of genes expression")
abline(h = 20e6, lty= 2)

# Get log2 counts per million
logcounts <- log2(countdata + 1)
