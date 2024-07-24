# install packages
install.packages(c("ggplot2", "data.table"))

# load packages
library(data.table)
library(ggplot2)

#set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Metagenomics/ATLAS/gene_catalogue/DESeq2/02_results/01_bf_nbf/")

# load file
data <- read.csv("DESeq2_BF_NBF.csv", sep = ",", row.names = 1)
data[["id"]] <- rownames(data)
head(data)

theme_minimal() |> theme_set()

# being more stringent here because there are a lot of them
sig <- data[data$padj < 0.01, ]

ggplot(sig) +
  aes(y=log2FoldChange, ymin=log2FoldChange - lfcSE, 
      ymax=log2FoldChange + lfcSE, x=reorder(id, -log2FoldChange), # can also do alphabetical, by p-value, etc.
      color=log2FoldChange<0) +                                         # could also color by baseMean
  geom_hline(yintercept=0, linetype="dashed") +
  geom_linerange() +
  geom_point() +
  labs(x="KO term", y="log2 fold-change") +
  guides(color=FALSE) +
  coord_flip()

# with different order
ggplot(sig) +
  aes(y=log2FoldChange, ymin=log2FoldChange - lfcSE, 
      ymax=log2FoldChange + lfcSE, x=reorder(id, -order), # can also do alphabetical, by p-value, etc.
      color=log2FoldChange<0) +                                         # could also color by baseMean
  geom_hline(yintercept=0, linetype="dashed") +
  geom_linerange() +
  geom_point() +
  labs(x="KO term", y="log2 fold-change") +
  guides(color=FALSE) +
  coord_flip()
