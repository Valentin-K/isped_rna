library(DESeq2)
library(tidyverse)
library("pasilla")


#---------- Extracting data from pasilla experiment
#-----------------------------------------------------


pasCts = system.file("extdata",
                     "pasilla_gene_counts.tsv",
                     package = "pasilla", mustWork = T)
pasAnno = system.file("extdata",
                      "pasilla_sample_annotation.csv",
                      package ="pasilla", mustWork = T)
cts = as.matrix(read.csv(pasCts, sep= "\t", row.names ="gene_id"))
coldata = read.csv (pasAnno, row.names=1)
coldata = coldata[,c("condition","type")]