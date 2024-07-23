# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)
library(microbiomeMarker)


######################################################################
##################### oral ###########################################
######################################################################

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/13_Permanova/")


#import data

oral <- import_qiime2(
  otu_qza = "../11_alpha_diversity/01_input/data_table_inf_dev_oral_bacteria.qza", 
  taxa_qza = "../11_alpha_diversity/01_input/taxonomy_inf_dev_oral_bacteria.qza",
  sam_tab = "01_input/metadata_inf_dev_oral_bacteria.txt")

detach("package:microbiomeMarker", unload = TRUE)  #unload microbiome maker again as it interferes with other commands later on


############################### birthmode ###########################################


# Filter the data based on the condition
pseq <- subset_samples(oral, timepoint == "M12")

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# plot
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "birthmode", size = 3)
print(p)

# samples x species as input
detach("package:microbiome", unload = TRUE) 
library(vegan)
permanova <- adonis2(t(otu) ~ birthmode,
                    data = meta, permutations=999, method = "bray")
permanova
write.csv(permanova, file = "02_output/permanova_oral_M12_birthmode.csv")




################# feeding ############################


# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)
library(microbiomeMarker)

# set working directory

setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/13_Permanova/")


#import data

oral <- import_qiime2(
  otu_qza = "../11_alpha_diversity/01_input/data_table_inf_dev_oral_bacteria.qza", 
  taxa_qza = "../11_alpha_diversity/01_input/taxonomy_inf_dev_oral_bacteria.qza",
  sam_tab = "01_input/metadata_inf_dev_oral_bacteria.txt")

detach("package:microbiomeMarker", unload = TRUE)  #unload microbiome maker again as it interferes with other commands later on


# Filter the data based on the condition
pseq <- subset_samples(oral, timepoint == "M12")
pseq <- subset_samples(pseq, breastfeeding =="x" | breastfeeding == "n")

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# plot
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "breastfeeding", size = 3)
print(p)

# samples x species as input
detach("package:microbiome", unload = TRUE) 
library(vegan)
permanova <- adonis2(t(otu) ~ breastfeeding,
                     data = meta, permutations=999, method = "bray")
permanova
write.csv(permanova, file = "02_output/permanova_oral_M12_feeding.csv")




################# timepoints per feeding-group ############################


# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)
library(microbiomeMarker)

# set working directory

setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/13_Permanova/")


#import data

oral <- import_qiime2(
  otu_qza = "../11_alpha_diversity/01_input/data_table_inf_dev_oral_bacteria.qza", 
  taxa_qza = "../11_alpha_diversity/01_input/taxonomy_inf_dev_oral_bacteria.qza",
  sam_tab = "01_input/metadata_inf_dev_oral_bacteria.txt")

detach("package:microbiomeMarker", unload = TRUE)  #unload microbiome maker again as it interferes with other commands later on

library(microbiome)
# Filter the data based on the condition
pseq <- subset_samples(oral, breastfeeding =="n")
pseq <- subset_samples(pseq, timepoint == "M01" |  timepoint =="M02")


# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# plot
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "timepoint", size = 3)
print(p)

# samples x species as input
detach("package:microbiome", unload = TRUE) 
library(vegan)
permanova <- adonis2(t(otu) ~ timepoint,
                     data = meta, permutations=999, method = "bray")
permanova
write.csv(permanova, file = "02_output/permanova_oral_breastfeeding_BF_M01_M02.csv")





######################################################################
##################### stool ##########################################
######################################################################

# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)
library(microbiomeMarker)

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/02_infants_development_stool_bacteria/13_Permanova/")


#import data

stool <- import_qiime2(
  otu_qza = "../10_alpha_diversity/01_input/data_table_inf_dev_stool_bacteria.qza", 
  taxa_qza = "../10_alpha_diversity/01_input/taxonomy_inf_dev_stool_bacteria.qza",
  sam_tab = "01_input/metadata_inf_dev_stool_bacteria.txt")

detach("package:microbiomeMarker", unload = TRUE)  #unload microbiome maker again as it interferes with other commands later on


############################### birthmode ###########################################


# Filter the data based on the condition
pseq <- subset_samples(stool, timepoint == "M12")

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# plot
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "birthmode", size = 3)
print(p)

# samples x species as input
detach("package:microbiome", unload = TRUE) 
library(vegan)
permanova <- adonis2(t(otu) ~ birthmode,
                     data = meta, permutations=999, method = "bray")
permanova
write.csv(permanova, file = "02_output/permanova_oral_M12_birthmode.csv")


################# feeding ############################

# Load libraries
library(microbiome)

# Filter the data based on the condition
pseq <- subset_samples(stool, timepoint == "M12")
pseq <- subset_samples(pseq, breastfeeding =="x" | breastfeeding == "n")

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# plot
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "breastfeeding", size = 3)
print(p)

# samples x species as input
detach("package:microbiome", unload = TRUE) 
library(vegan)
permanova <- adonis2(t(otu) ~ breastfeeding,
                     data = meta, permutations=999, method = "bray")
permanova
write.csv(permanova, file = "02_output/permanova_stool_M12_feeding.csv")



################# timepoints per feeding ############################

# Load libraries
library(microbiome)

# Filter the data based on the condition
pseq <- subset_samples(stool, breastfeeding =="n")
pseq <- subset_samples(pseq, timepoint == "M12" | timepoint =="M10")


# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# plot
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "timepoint", size = 3)
print(p)

# samples x species as input
detach("package:microbiome", unload = TRUE) 
library(vegan)
permanova <- adonis2(t(otu) ~ timepoint,
                     data = meta, permutations=999, method = "bray")
permanova
write.csv(permanova, file = "02_output/permanova_stool_NBF_M12_M10.csv")
