# install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")

# load libraries
library(phyloseq)
library(ggplot2)
library(decontam)

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/infant_development/Amplikon_Analysis/04_infants_development_stool_archaea/06_decontam")
# load OTU file 
# delete first line, delete # in #OTU ID
# deleted NC that do not fit to stool
# Neg_Swab1	NegA_Swab3 
data_table <- read.delim("infants_development_stool_archaea_feature_table_from_biom_4decontam.txt")
View(data_table)

# load metadata file
metadata <- read.delim2("metadata_infants_development_stool_archaea_4decontam.txt")
View(metadata)

#create feature_table_universal
data_table_universal <- read.table("infants_development_stool_archaea_feature_table_from_biom_4decontam.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

#create map_samples_universal
map_samples_universal <- read.table("metadata_infants_development_stool_archaea_4decontam.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

#create ps_d_universal with phyloseq
ps_d_universal <- phyloseq(otu_table(data_table_universal, taxa_are_rows = TRUE), sample_data(map_samples_universal))
ps_d_universal

# output
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2846 taxa and 461 samples ]
# sample_data() Sample Data:       [ 461 samples by 1 sample variables ]

sample_data(ps_d_universal)$is.neg <- sample_data(ps_d_universal)$taxon == "control"
contamdf.prev <- isContaminant(ps_d_universal, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)

# output
# FALSE  TRUE 
# 2788    58

# not needed
# head(which(contamdf.prev$contaminant))


# write contaminant ASVs in table
write.table(which(contamdf.prev$contaminant), file="infants_development_stool_archaea_iscontaminant_prevalence_05.tsv")

#create new "clean" feature table without contaminants
data_table_clean_universal <- data_table_universal[!contamdf.prev$contaminant, ]
write.table(file = "infants_development_stool_archaea_iscontaminant_prevalence_05_clean.tsv", x = data.frame(data_table_clean_universal), sep = "\t", row.names = TRUE, col.names = TRUE)



### try different thresholds -> 0.3
contamdf.prev <- isContaminant(ps_d_universal, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev$contaminant)

# output
# FALSE  TRUE 
# 2799    47

# not needed
# head(which(contamdf.prev$contaminant))


# write contaminant ASVs in table
write.table(which(contamdf.prev$contaminant), file="infants_development_stool_archaea_iscontaminant_prevalence_03.tsv")

#create new "clean" feature table without contaminants
data_table_clean_universal <- data_table_universal[!contamdf.prev$contaminant, ]
write.table(file = "infants_development_stool_archaea_iscontaminant_prevalence_03_clean.tsv", x = data.frame(data_table_clean_universal), sep = "\t", row.names = TRUE, col.names = TRUE)



### try different thresholds  -> 0.1
contamdf.prev <- isContaminant(ps_d_universal, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev$contaminant)

# output
# FALSE  TRUE 
# 2812    34

# not needed
# head(which(contamdf.prev$contaminant))


# write contaminant ASVs in table
write.table(which(contamdf.prev$contaminant), file="infants_development_stool_archaea_iscontaminant_prevalence_01.tsv")

#create new "clean" feature table without contaminants
data_table_clean_universal <- data_table_universal[!contamdf.prev$contaminant, ]
write.table(file = "infants_development_stool_archaea_iscontaminant_prevalence_01_clean.tsv", x = data.frame(data_table_clean_universal), sep = "\t", row.names = TRUE, col.names = TRUE)
