# DE Analysis of RNA-Seq data

#Import cts and coldata
#cts countains the count matrix
#coldata contains info about the genotype and treatment for each sample
cts <- read.table(file = 'cts.tsv', sep = '\t', header = TRUE)
coldata <- read.table(file = 'coldata.txt', sep = '\t', header = TRUE)

#Ens_gene names as row_names
rownames(cts)=cts[1:38293,1]
cts <- cts[1:38293,-1]
#Samples as row_names
rownames(coldata)=coldata[1:24,1]
coldata <- coldata[1:24,-1]

###WT_TRT VS WT_CTR
cts1 <- cts[,c(8,9,12,19,20,21, 2,3,6,13,14,15)]
coldata1 <-coldata[c(8,9,12,19,20,21, 2,3,6,13,14,15),]

###KO_TRT VS KO_CTR
cts2 <- cts[,c(7,10,11,22,23,24, 1,4,5,16,17,18)]
coldata2 <- coldata[c(7,10,11,22,23,24, 1,4,5,16,17,18),]

###KO_CTR VS WT_CTR  
cts3 <- cts[,c(1,4,5,16,17,18, 2,3,6,13,14,15)]
coldata3 <- coldata[c(1,4,5,16,17,18, 2,3,6,13,14,15),]

###KO_TRT VS WT_TRT  
cts4 <- cts[,c(7,10,11,22,23,24, 8,9,12,19,20,21)]
coldata4 <- coldata[c(7,10,11,22,23,24, 8,9,12,19,20,21),]





#Delete genes with counts for all samples < 12
#rowSums(...matrix...)
cts1 <- cts1[(rowSums(cts1)>12),]
cts2 <- cts2[(rowSums(cts2)>12),]
cts3 <- cts3[(rowSums(cts3)>12),]
cts4 <- cts4[(rowSums(cts4)>12),]



#For the DE analysis we use DESeq2
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)


dds1 <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = coldata1,
                              design = ~ Treatment)
dds1
dds1 <- DESeq(dds1)
res1 <- results(dds1)
write.table(res1, sep = "\t", file = "WT_TRT VS WT_CTR")



dds2 <- DESeqDataSetFromMatrix(countData = cts2,
                               colData = coldata2,
                               design = ~ Treatment)
dds2
dds2 <- DESeq(dds2)
res2 <- results(dds2)
write.table(res2, sep = "\t", file = "KO_TRT VS KO_CTR")



dds3 <- DESeqDataSetFromMatrix(countData = cts3,
                               colData = coldata3,
                               design = ~ Genotype)
dds3
dds3 <- DESeq(dds3)
res3 <- results(dds3)
write.table(res3, sep = "\t", file = "KO_CTR VS WT_CTR")



dds4 <- DESeqDataSetFromMatrix(countData = cts4,
                               colData = coldata4,
                               design = ~ Genotype)
dds4
dds4 <- DESeq(dds4)
res4 <- results(dds4)
write.table(res4, sep = "\t", file = "KO_TRT VS WT_TRT")
