library(DESeq2)
library(tidyverse)
counts <- read.delim("counts.tsv", stringsAsFactors=F)

countdata <- counts
keep <- rowSums(countdata) >5
countdata <- countdata[keep,]
dim(countdata)

#Bar plot of library sizes
librarySizes <- colSums(countdata)
barplot (libraarySizes,
         names = names(librarySizes),
         las=2,
         main="Barplot of library sizes")
abline(h = 20e6, lty= 2)
