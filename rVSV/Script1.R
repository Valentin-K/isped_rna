library(DESeq2)
library(tidyverse)
counts <- read.delim("counts.tsv", stringsAsFactors=F)

countdata <- counts %>%
  column_to_rownames("X") %>% # turn the geneid column into rownames
  as.matrix()
head(countdata)


keep <- rowSums(countdata, na.rm=F) >5
countdata <- countdata[keep,]
dim(countdata)

#Bar plot of library sizes
librarySizes <- colSums(countdata)
barplot (librarySizes,
         names = names(librarySizes),
         las=2,
         main="Barplot of library sizes")
abline(h = 20e6, lty= 2)
