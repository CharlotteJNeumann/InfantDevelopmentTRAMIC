#How to analyze your gene catalog output from ATLAS 
#Normalization output exporting
#Input here is the unnormalized data from Silas script
# Method used for normalization in Deseq: relative log expression (RLE)

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Metagenomics/ATLAS/gene_catalogue/DESeq2/")

#Install Deseq2 first
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#######################################################
########## birth mode #################################
#######################################################

#Read the csv file (gene_annotation) and import as as matrix
gene_data <- as.matrix(read.csv("01_input/annotation_genecounts_unnormalized_M01.csv", header = TRUE, sep = ",", dec = ".", row.names="gene_id"))
View(gene_data)

#Import the metadata file 

metadata <- read.csv("01_input/metadata_M01.csv", header = TRUE, sep = ",", dec = ".", row.names="Sample_ID")
View(metadata)

#Define the columns and metadata (condition or type or age, or study, etc)
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

metadata$birthmode <- factor(metadata$birthmode)

levels(metadata$birthmode)


#With the count matrix, gene_data, and the sample information, metadata, we can construct a DESeqDataSet:
#count data is the gene_data (gene count)
#colData is the metadata file
#design is the grouping
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = gene_data,
                              colData = metadata,
                              design = ~ birthmode)
dds

#To get the normalized counts use one of the followings. Both outputs will be the same. 
# 1)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="02_results/02_sp_CS/normalized_counts_M12.csv")

# 2)
#dds <- DESeq(dds)
#sizeFactors(dds)
#normalized_counts_DESeq <- counts(dds, normalized=TRUE)
#write.csv(normalized_counts_DESeq, file="normalized_counts.csv")

#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)

#Pre-filtering
#present in at least 3 samples (A recommendation for the minimal number of samples is to specify the smallest group size: in my case number of high methane producers
#The count of 10 is a reasonable choice for bulk RNA-seq (min number of reads) #number of methane producers in adults=32 
smallestGroupSize <- 2 
keep <- rowSums(counts(dds) >= 2) >= smallestGroupSize
dds <- dds[keep,]
#choose a reference level 
# dds$birthmode <- relevel(dds$birthmode, ref = "sp")
 dds$birthmode <- relevel(dds$birthmode, ref = "sp")
#------------------------------------------------------------------

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

#Exploring and exporting results
plotMA(res, ylim=c(-2,2))


#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]

#Exporting results to CSV files, here you will see the diff abundant genes
write.csv(as.data.frame(resOrdered), file="02_results/02_sp_CS/Deseq2_M12_results.csv")





#######################################################
############# feeding #################################
#######################################################


#Read the csv file (gene_annotation) and import as as matrix
gene_data <- as.matrix(read.csv("01_input/annotation_genecounts_unnormalized_M12.csv", header = TRUE, sep = ",", dec = ".", row.names="gene_id"))
View(gene_data)

#Import the metadata file 

metadata <- read.csv("01_input/metadata_M12.csv", header = TRUE, sep = ",", dec = ".", row.names="Sample_ID")
View(metadata)
# filter NAs
#not needed for M12
metadata <- subset(metadata, feeding =="BF" | feeding =="NBF")
View(metadata)
#Define the columns and metadata (condition or type or age, or study, etc)
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

metadata$feeding <- factor(metadata$feeding)
levels(metadata$feeding)



#With the count matrix, gene_data, and the sample information, metadata, we can construct a DESeqDataSet:
#count data is the gene_data (gene count)
#colData is the metadata file
#design is the grouping
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = gene_data,
                              colData = metadata,
                              design = ~ feeding)
dds

#To get the normalized counts use one of the followings. Both outputs will be the same. 
# 1)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="02_results/01_bf_nbf/normalized_counts_M12.csv")

# 2)
#dds <- DESeq(dds)
#sizeFactors(dds)
#normalized_counts_DESeq <- counts(dds, normalized=TRUE)
#write.csv(normalized_counts_DESeq, file="normalized_counts.csv")

#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)

#Pre-filtering
#present in at least 3 samples (A recommendation for the minimal number of samples is to specify the smallest group size: in my case number of high methane producers
#The count of 10 is a reasonable choice for bulk RNA-seq (min number of reads) #number of methane producers in adults=32 
smallestGroupSize <- 2 
keep <- rowSums(counts(dds) >= 2) >= smallestGroupSize
dds <- dds[keep,]
#choose a reference level 
dds$feeding <- relevel(dds$feeding, ref = "BF")
#------------------------------------------------------------------

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

#Exploring and exporting results
plotMA(res, ylim=c(-2,2))


#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]

#Exporting results to CSV files, here you will see the diff abundant genes
write.csv(as.data.frame(resOrdered), file="02_results/01_bf_nbf/Deseq2_M12_results.csv")



#######################################################
########## time point #################################
#######################################################

#Read the csv file (gene_annotation) and import as as matrix
gene_data <- as.matrix(read.csv("01_input/annotation_genecounts_unnormalized_M01_M06_M12.csv", header = TRUE, sep = ",", dec = ".", row.names="gene_id"))
View(gene_data)

#Import the metadata file 

metadata <- read.csv("01_input/metadata.csv", header = TRUE, sep = ",", dec = ".", row.names="Sample_ID")
View(metadata)

#Define the columns and metadata (condition or type or age, or study, etc)
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

metadata$tp <- factor(metadata$tp)

levels(metadata$tp)


#With the count matrix, gene_data, and the sample information, metadata, we can construct a DESeqDataSet:
#count data is the gene_data (gene count)
#colData is the metadata file
#design is the grouping
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = gene_data,
                              colData = metadata,
                              design = ~ tp)
dds

#To get the normalized counts use one of the followings. Both outputs will be the same. 
# 1)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="02_results/03_tps/normalized_counts_M01_M06_M12.csv")

# 2)
#dds <- DESeq(dds)
#sizeFactors(dds)
#normalized_counts_DESeq <- counts(dds, normalized=TRUE)
#write.csv(normalized_counts_DESeq, file="normalized_counts.csv")

#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)

#Pre-filtering
#present in at least 3 samples (A recommendation for the minimal number of samples is to specify the smallest group size: in my case number of high methane producers
#The count of 10 is a reasonable choice for bulk RNA-seq (min number of reads) #number of methane producers in adults=32 
smallestGroupSize <- 2 
keep <- rowSums(counts(dds) >= 2) >= smallestGroupSize
dds <- dds[keep,]
#choose a reference level 
# dds$birthmode <- relevel(dds$birthmode, ref = "sp")
dds$tp <- relevel(dds$tp, ref = "M01", "M06") #change for M06 M12 and M01 M12
?relevel
#------------------------------------------------------------------

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

#Exploring and exporting results
plotMA(res, ylim=c(-2,2))


#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]

#Exporting results to CSV files, here you will see the diff abundant genes
write.csv(as.data.frame(resOrdered), file="02_results/03_tps/Deseq2_M01_M06_results_referenceM01.csv")




#######################################################
########## time point per feeding type ################
#######################################################

#Read the csv file (gene_annotation) and import as as matrix
gene_data <- as.matrix(read.csv("01_input/annotation_genecounts_unnormalized_M01_M06_M12_NBF.csv", header = TRUE, sep = ",", dec = ".", row.names="gene_id"))
View(gene_data)

#Import the metadata file 

metadata <- read.csv("01_input/metadata.csv", header = TRUE, sep = ",", dec = ".", row.names="Sample_ID")
View(metadata)

# filter feeding
metadata <- subset(metadata, feeding =="NBF")
View(metadata)
#Define the columns and metadata (condition or type or age, or study, etc)
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

metadata$feeding <- factor(metadata$feeding)
levels(metadata$feeding)

#Define the columns and metadata (condition or type or age, or study, etc)
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

metadata$tp <- factor(metadata$tp)

levels(metadata$tp)


#With the count matrix, gene_data, and the sample information, metadata, we can construct a DESeqDataSet:
#count data is the gene_data (gene count)
#colData is the metadata file
#design is the grouping
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = gene_data,
                              colData = metadata,
                              design = ~ tp)
dds

#To get the normalized counts use one of the followings. Both outputs will be the same. 
# 1)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="02_results/04_tps_BF_NBF/normalized_counts_M01_M06_M12_NBF.csv")

# 2)
#dds <- DESeq(dds)
#sizeFactors(dds)
#normalized_counts_DESeq <- counts(dds, normalized=TRUE)
#write.csv(normalized_counts_DESeq, file="normalized_counts.csv")

#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)

#Pre-filtering
#present in at least 3 samples (A recommendation for the minimal number of samples is to specify the smallest group size: in my case number of high methane producers
#The count of 10 is a reasonable choice for bulk RNA-seq (min number of reads) #number of methane producers in adults=32 
smallestGroupSize <- 2 
keep <- rowSums(counts(dds) >= 2) >= smallestGroupSize
dds <- dds[keep,]
#choose a reference level 
dds$tp <- relevel(dds$tp, ref = "M12") #change for M06 M12 and M01 M12
levels(dds$tp)
?relevel
#------------------------------------------------------------------

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

#Exploring and exporting results
plotMA(res, ylim=c(-2,2))


#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]

#Exporting results to CSV files, here you will see the diff abundant genes
write.csv(as.data.frame(resOrdered), file="02_results/04_tps_BF_NBF/Deseq2_NBF_M06_M12_results_referenceM12.csv")


