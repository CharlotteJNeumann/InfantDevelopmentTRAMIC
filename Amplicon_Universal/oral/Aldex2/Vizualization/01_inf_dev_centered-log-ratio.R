###################################################
########## CLR transformation #####################
###################################################

# I usually run Aldex2 in qiime on the cluster
# for vizualizing my data in R in boxplots I need CLR transformed data
# so this script performes CLR (centered log ration) transformation
# as input: genus collapsed data table
# output used for boxplots


# load packages
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(purrr)
library(broom)
library(pheatmap)
library(plotly)
library(microbiome)
library(ggbeeswarm)
library(knitr)
library(ALDEx2)

########################################################################
########################### ORAL #######################################
########################################################################

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/12_Aldex2/03_vizualization")

#read in file
abundance_file_TRAMIC = "data_table_infants_development_oral_bacteria_genus.txt"

#For relative abundance
D_TRAMIC <- read_tsv(abundance_file_TRAMIC, show_col_types = FALSE) %>%
  column_to_rownames(var = "OTU") %>%
  as.data.frame()

# calculate relative abundance
rel_ab_TRAMIC <- sweep(D_TRAMIC, 1, colSums(D_TRAMIC),`/`)

# transform counts with centered log ratio
clr_data_TRAMIC <- transform(rel_ab_TRAMIC, transform = "clr")

write.csv(clr_data_TRAMIC, file ="clr_data_collapsed_genus.csv")






########################################################################
########################### STOOL ######################################
########################################################################

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/02_infants_development_stool_bacteria/12_Aldex2/03_vizualization")

#read in file
abundance_file_TRAMIC = "datafile_infants_development_stool_bacteria_genus.txt"

#For relative abundance
D_TRAMIC <- read_tsv(abundance_file_TRAMIC, show_col_types = FALSE) %>%
  column_to_rownames(var = "...1") %>%
  as.data.frame()

# calculate relative abundance
rel_ab_TRAMIC <- sweep(D_TRAMIC, 1, colSums(D_TRAMIC),`/`)

# transform counts with centered log ratio
clr_data_TRAMIC <- transform(rel_ab_TRAMIC, transform = "clr")

write.csv(clr_data_TRAMIC, file ="clr_data_genus.csv")
