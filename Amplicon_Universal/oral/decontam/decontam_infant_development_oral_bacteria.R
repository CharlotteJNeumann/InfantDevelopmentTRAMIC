# install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")

# load libraries
library(phyloseq)
library(ggplot2)
library(decontam)

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/06_decontam")

# load OTU file
# v2: removed negative controls which do not fit to oral:
# B_NC_iS1, B_NC_iS2, B_NC_iS3, B_NC_i01_S4 B_Neg_stool1	B_Neg_Stool1	B_Neg_stool2, B_NegStool3
data_table <- read.delim("infants_development_oral_bacteria_feature-table_from_biom_4decontam_v2.txt")
View(data_table)

# load metadata file
metadata <- read.delim2("Metadata_infants_development_oral_bacteria_4decontam.txt")
View(metadata)

#create feature_table_universal
data_table_universal <- read.table("infants_development_oral_bacteria_feature-table_from_biom_4decontam_v2.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

#create map_samples_universal
map_samples_universal <- read.table("Metadata_infants_development_oral_bacteria_4decontam.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

#create ps_d_universal with phyloseq
ps_d_universal <- phyloseq(otu_table(data_table_universal, taxa_are_rows = TRUE), sample_data(map_samples_universal))
ps_d_universal

# output
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 11821 taxa and 430 samples ]
# sample_data() Sample Data:       [ 430 samples by 1 sample variables ]

sample_data(ps_d_universal)$is.neg <- sample_data(ps_d_universal)$taxon == "control"
contamdf.prev <- isContaminant(ps_d_universal, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)

# output
# FALSE  TRUE 
# 11511   310 

head(which(contamdf.prev$contaminant))
# output
# [1] 21 29 32 35 49 51

# write contaminant ASVs in table
write.table(which(contamdf.prev$contaminant), file="infants_development_oral_bacteria_v2_iscontaminant_prevalence_05.tsv")

#create new "clean" feature table without contaminants
data_table_clean_universal <- data_table_universal[!contamdf.prev$contaminant, ]
write.table(file = "infants_development_oral_bacteria_v2_iscontaminant_prevalence_05_clean.tsv", x = data.frame(data_table_clean_universal), sep = "\t", row.names = TRUE, col.names = TRUE)


### try different thresholds
contamdf.prev <- isContaminant(ps_d_universal, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev$contaminant)

# output
# FALSE  TRUE 
# 11632   189  

head(which(contamdf.prev$contaminant))
# output
# [1] 29 35 70 71 78 81

# write contaminant ASVs in table
write.table(which(contamdf.prev$contaminant), file="infants_development_oral_bacteria_v2_iscontaminant_prevalence_01.tsv")

#create new "clean" feature table without contaminants
data_table_clean_universal <- data_table_universal[!contamdf.prev$contaminant, ]
write.table(file = "infants_development_oral_bacteria_v2_iscontaminant_prevalence_01_clean.tsv", x = data.frame(data_table_clean_universal), sep = "\t", row.names = TRUE, col.names = TRUE)


### try different thresholds
contamdf.prev <- isContaminant(ps_d_universal, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev$contaminant)

# output
# FALSE  TRUE 
# 11540   281  

head(which(contamdf.prev$contaminant))
# output
# [1] 29 35 59 67 70 71

# write contaminant ASVs in table
write.table(which(contamdf.prev$contaminant), file="infants_development_oral_bacteria_v2_iscontaminant_prevalence_03.tsv")

#create new "clean" feature table without contaminants
data_table_clean_universal <- data_table_universal[!contamdf.prev$contaminant, ]
write.table(file = "infants_development_oral_bacteria_v2_iscontaminant_prevalence_03_clean.tsv", x = data.frame(data_table_clean_universal), sep = "\t", row.names = TRUE, col.names = TRUE)
