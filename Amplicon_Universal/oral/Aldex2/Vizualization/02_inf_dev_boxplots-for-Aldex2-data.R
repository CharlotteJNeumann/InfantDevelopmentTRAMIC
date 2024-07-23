#################################################################################

####################### R script BOX PLOTS of ALDEX2 ################################

#################################################################################

#################################################################################
################################ ORAL ###########################################
##################################################################################



# load library
library(ggplot2)
library(dplyr)
library(reshape2)
library(phyloseq)

# use clr transofrmed data table (previous script) as input

setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/12_Aldex2/03_vizualization")

#read in metadata
metadata = read.table("Metadata_infants_development_oral_bacteria.txt",row.names = 1)
#read in otu table
otu_table = read.csv("clr_data_genus_Aldex2_sgnf_diff.csv", sep = ",", row.names = 1)

# Merge metadata with OTU table
merged_data <- merge(metadata, otu_table, by = "row.names")

# Reshape data into long format
merged_data_long <- melt(merged_data, id.vars = c("Row.names", colnames(metadata)))

#rename the variables
names(merged_data_long)[names(merged_data_long) == "Row.names"] <- "sample"
names(merged_data_long)[names(merged_data_long) == "V2"] <- "individual"
names(merged_data_long)[names(merged_data_long) == "V3"] <- "timepoint"
names(merged_data_long)[names(merged_data_long) == "V4"] <- "feeding"
names(merged_data_long)[names(merged_data_long) == "variable"] <- "genus"

## create subsets for time and genera that are significantly differentially abundant in Aldex2

# M02
M02 <- subset(merged_data_long, timepoint == "M02")
M02 <- subset(M02, feeding == "n"| feeding == "x")
M02_sgnf <- subset(M02, genus == "Actinomyces" |genus == "Prevotella" | genus == "Granulicatella" | 
                     genus =="Lachnoanaerobaculum")

# M03
M03 <- subset(merged_data_long, timepoint == "M03")
M03 <- subset(M03, feeding == "n"| feeding == "x")
M03_sgnf <- subset(M03, genus == "Actinomyces" |genus == "Prevotella" | genus == "Granulicatella" | 
                     genus =="Veillonella" |genus =="Streptococcus" | genus =="Actinomyces" | 
                     genus =="Lachnoanaerobaculum")

#M04
M04 <- subset(merged_data_long, timepoint == "M04")
M04 <- subset(M04, feeding == "n"| feeding == "x")
M04_sgnf <- subset(M04, genus == "Prevotella" | genus =="Lachnoanaerobaculum" | genus == "Campylobacter" | 
                     genus =="Leptotrichia" | genus =="Veillonella" | genus =="Granulicatella" | 
                     genus =="Streptococcus" | genus == "Actinomyces" |genus =="Atopobium" |
                     genus == "Neisseria" | genus == "Haemophilus" | genus == "Fusobacterium")

#M05
M05 <- subset(merged_data_long, timepoint == "M05")
M05 <- subset(M05, feeding == "n"| feeding == "x")
M05_sgnf <- subset(M05, genus == "Fusobacterium" | genus =="Streptococcus" | genus == "Solobacterium" | 
                     genus =="Leptotrichia" | genus =="Granulicatella" | genus =="Prevotella" | 
                     genus == "Lachnoanaerobaculum" | genus == "Neisseria")

#M06
M06 <- subset(merged_data_long, timepoint == "M06")
M06 <- subset(M06, feeding == "n"| feeding == "x")
M06_sgnf <- subset(M06, genus == "Lachnoanaerobaculum" | genus =="Leptotrichia" | genus == "Streptococcus" | 
                     genus =="Neisseria" | genus =="Prevotella")

#M07
M07 <- subset(merged_data_long, timepoint == "M07")
M07 <- subset(M07, feeding == "n"| feeding == "x")
M07_sgnf <- subset(M07, genus == "Lachnoanaerobaculum" | genus =="Neisseria" | genus =="Fusobacterium" | 
                     genus == "unclassified_Pasteurellaceae" | genus =="Prevotella")

#M08
M08 <- subset(merged_data_long, timepoint == "M08")
M08 <- subset(M08, feeding == "n"| feeding == "x")
M08_sgnf <- subset(M08, genus =="Fusobacterium" | genus =="Prevotella" | genus == "Gemella" |
                     genus == "Solobacterium")

#M09
M09 <- subset(merged_data_long, timepoint == "M09")
M09 <- subset(M09, feeding == "n"| feeding == "x")
M09_sgnf <- subset(M09, genus == "Veillonella" | genus =="Streptococcus" | genus =="Neisseria")

#M10
M10 <- subset(merged_data_long, timepoint == "M10")
M10 <- subset(M10, feeding == "n"| feeding == "x")
M10_sgnf <- subset(M10, genus == "Solobacterium" | genus =="Gemella" | genus =="Veillonella")

#M11
M11 <- subset(merged_data_long, timepoint == "M11")
M11 <- subset(M11, feeding == "n"| feeding == "x")
M11_sgnf <- subset(M11, genus =="Neisseria")


### Plot, print and save using ggplot2
## M02
#plot
plot_M02_sgnf <- ggplot(M02_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M02_sgnf)

#save the plot
ggsave("M02_sgnf.png", plot_M02_sgnf, width = 15, height = 15, units = "cm")
ggsave("M02_sgnf.svg", plot_M02_sgnf, width = 15, height = 15, units = "cm")

### Plot, print and save using ggplot2
## M03
#plot
plot_M03_sgnf <- ggplot(M03_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M03_sgnf)

#save the plot
ggsave("M03_sgnf.png", plot_M03_sgnf, width = 20, height = 15, units = "cm")
ggsave("M03_sgnf.svg", plot_M03_sgnf, width = 20, height = 15, units = "cm")

### Plot, print and save using ggplot2
## M04
#plot
plot_M04_sgnf <- ggplot(M04_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M04_sgnf)

#save the plot
ggsave("M04_sgnf.png", plot_M04_sgnf, width = 35, height = 15, units = "cm")
ggsave("M04_sgnf.svg", plot_M04_sgnf, width = 35, height = 15, units = "cm")

### Plot, print and save using ggplot2
## M05
#plot
plot_M05_sgnf <- ggplot(M05_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M05_sgnf)

#save the plot
ggsave("M05_sgnf.png", plot_M05_sgnf, width = 25, height = 15, units = "cm")
ggsave("M05_sgnf.svg", plot_M05_sgnf, width = 25, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M06
#plot
plot_M06_sgnf <- ggplot(M06_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M06_sgnf)

#save the plot
ggsave("M06_sgnf.png", plot_M06_sgnf, width = 17.5, height = 15, units = "cm")
ggsave("M06_sgnf.svg", plot_M06_sgnf, width = 17.5, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M07
#plot
plot_M07_sgnf <- ggplot(M07_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M07_sgnf)

#save the plot
ggsave("M07_sgnf.png", plot_M07_sgnf, width = 15, height = 15, units = "cm")
ggsave("M07_sgnf.svg", plot_M07_sgnf, width = 15, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M08
#plot
plot_M08_sgnf <- ggplot(M08_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M08_sgnf)

#save the plot
ggsave("M08_sgnf.png", plot_M08_sgnf, width = 15, height = 15, units = "cm")
ggsave("M08_sgnf.svg", plot_M08_sgnf, width = 15, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M09
#plot
plot_M09_sgnf <- ggplot(M09_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M09_sgnf)

#save the plot
# genus: 3
ggsave("M09_sgnf.png", plot_M09_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M09_sgnf.svg", plot_M09_sgnf, width = 12.5, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M10
#plot
plot_M10_sgnf <- ggplot(M10_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M10_sgnf)

#save the plot
ggsave("M10_sgnf.png", plot_M10_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M10_sgnf.svg", plot_M10_sgnf, width = 12.5, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M11
#plot
plot_M11_sgnf <- ggplot(M11_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M11_sgnf)

#save the plot
ggsave("M11_sgnf.png", plot_M11_sgnf, width = 7.5, height = 15, units = "cm")
ggsave("M11_sgnf.svg", plot_M11_sgnf, width = 7.5, height = 15, units = "cm")



############################# p-values ###################
# NOT WORKING!!!

### Plot, print and save using ggplot2
## M03

p_values <- matrix(c(0.05, 0.01, 0.001), nrow = 3, ncol = 3,
                   dimnames = list(c("A vs B", "A vs C", "B vs C"),
                                   c("p_value")))
# Define the positions where significance indicators will be placed (example)
signif_positions_M03 <- c(1.5, 2.5, 3.5, 0.001, 0.04, 0.00004)  # Adjust as needed based on your plot
#plot
plot_M03_sgnf <- ggplot(M03_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()
# Add significance indicators (*)
geom_text(data = data.frame(x = signif_positions_M03, y = max(M03_sgnf$value), label = "*"),
          aes(x = x, y = y, label = label), 
          size = 6, vjust = -0.5, hjust = 0)

#print
print(plot_M03_sgnf)

#save the plot
ggsave("M03_sgnf.png", plot_M03_sgnf, width = 20, height = 15, units = "cm")
ggsave("M03_sgnf.svg", plot_M03_sgnf, width = 20, height = 15, units = "cm")









########################################################################################
################################# STOOL ################################################
#######################################################################################


# load library
library(ggplot2)
library(dplyr)
library(reshape2)
library(phyloseq)

# use clr transformed data table (previous script) as input

setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/02_infants_development_stool_bacteria/12_Aldex2/03_vizualization")

#read in metadata
metadata = read.table("Metadata_infants_development_stool_bacteria.txt", header = TRUE, row.names = 1)
#metadata = read.table("Metadata_infants_development_stool_bacteria.txt", row.names = 1)

#read in otu table (transposed!!!)
otu_table = read.csv("stool_clr_data_genus_Aldex2_sgnf_diff.csv", sep = ",", row.names = 1)

# Merge metadata with OTU table
merged_data <- merge(metadata, otu_table, by = "row.names")

# Reshape data into long format
merged_data_long <- melt(merged_data, id.vars = c("Row.names", colnames(metadata)))

#rename the variables
#names(merged_data_long)[names(merged_data_long) == "Row.names"] <- "sample"
#names(merged_data_long)[names(merged_data_long) == "V2"] <- "individual"
#names(merged_data_long)[names(merged_data_long) == "V3"] <- "timepoint"
#names(merged_data_long)[names(merged_data_long) == "V4"] <- "feeding"
names(merged_data_long)[names(merged_data_long) == "variable"] <- "genus"

## create subsets for time and genera that are significantly differentially abundant in Aldex2

# M01
M01 <- subset(merged_data_long, timepoint == "M01")
M01 <- subset(M01, feeding == "n"| feeding == "x")
M01_sgnf <- subset(M01, genus == "Enterococcus")

# M02
M02 <- subset(merged_data_long, timepoint == "M02")
M02 <- subset(M02, feeding == "n"| feeding == "x")
M02_sgnf <- subset(M02, genus == "Enterococcus" |genus == "Gemella" | genus == "Intestinibacter")

# M03
M03 <- subset(merged_data_long, timepoint == "M03")
M03 <- subset(M03, feeding == "n"| feeding == "x")
M03_sgnf <- subset(M03, genus == "Intestinibacter")

#M04
M04 <- subset(merged_data_long, timepoint == "M04")
M04 <- subset(M04, feeding == "n"| feeding == "x")
M04_sgnf <- subset(M04, genus == "Intestinibacter")

#M05
M05 <- subset(merged_data_long, timepoint == "M05")
M05 <- subset(M05, feeding == "n"| feeding == "x")
M05_sgnf <- subset(M05, genus == "Intestinibacter")

#M08
M08 <- subset(merged_data_long, timepoint == "M08")
M08 <- subset(M08, feeding == "n"| feeding == "x")
M08_sgnf <- subset(M08, genus =="Intestinibacter" | genus =="Lactobacillus" | genus == "Coprobacillus" )

#M09
M09 <- subset(merged_data_long, timepoint == "M09")
M09 <- subset(M09, feeding == "n"| feeding == "x")
M09_sgnf <- subset(M09, genus == "Lactobacillus" | genus == "Blautia" | genus == "Ruminococcus_gnavus_group")

#M10
M10 <- subset(merged_data_long, timepoint == "M10")
M10 <- subset(M10, feeding == "n"| feeding == "x")
M10_sgnf <- subset(M10, genus == "Lactobacillus" | genus =="EscherichiaShigella")

#M11
M11 <- subset(merged_data_long, timepoint == "M11")
M11 <- subset(M11, feeding == "n"| feeding == "x")
M11_sgnf <- subset(M11, genus =="Lactobacillus" | genus == "Blautia" | genus == "Megasphaera")

#M12
M12 <- subset(merged_data_long, timepoint == "M12")
M12 <- subset(M12, feeding == "n"| feeding == "x")
M12_sgnf <- subset(M12, genus =="Blautia" | genus == "Megasphaera" | genus == "Enterococcus")


### Plot, print and save using ggplot2

## M01
#plot
plot_M01_sgnf <- ggplot(M01_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M01_sgnf)

#save the plot
# genus: 1
ggsave("M01_sgnf.png", plot_M01_sgnf, width = 7.5, height = 15, units = "cm")
ggsave("M01_sgnf.svg", plot_M01_sgnf, width = 7.5, height = 15, units = "cm")


## M02
#plot
plot_M02_sgnf <- ggplot(M02_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M02_sgnf)

#save the plot
# genus: 3
ggsave("M02_sgnf.png", plot_M02_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M02_sgnf.svg", plot_M02_sgnf, width = 12.5, height = 15, units = "cm")

### Plot, print and save using ggplot2
## M03
#plot
plot_M03_sgnf <- ggplot(M03_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M03_sgnf)

#save the plot
# genus: 1
ggsave("M03_sgnf.png", plot_M03_sgnf, width = 7.5, height = 15, units = "cm")
ggsave("M03_sgnf.svg", plot_M03_sgnf, width = 7.5, height = 15, units = "cm")

### Plot, print and save using ggplot2
## M04
#plot
plot_M04_sgnf <- ggplot(M04_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M04_sgnf)

#save the plot
# genus: 1
ggsave("M04_sgnf.png", plot_M04_sgnf, width = 7.5, height = 15, units = "cm")
ggsave("M04_sgnf.svg", plot_M04_sgnf, width = 7.5, height = 15, units = "cm")

### Plot, print and save using ggplot2
## M05
#plot
plot_M05_sgnf <- ggplot(M05_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M05_sgnf)

#save the plot
# genus: 1
ggsave("M05_sgnf.png", plot_M05_sgnf, width = 7.5, height = 15, units = "cm")
ggsave("M05_sgnf.svg", plot_M05_sgnf, width = 7.5, height = 15, units = "cm")



### Plot, print and save using ggplot2
## M08
#plot
plot_M08_sgnf <- ggplot(M08_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M08_sgnf)

#save the plot
# genus: 3
ggsave("M08_sgnf.png", plot_M08_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M08_sgnf.svg", plot_M08_sgnf, width = 12.5, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M09
#plot
plot_M09_sgnf <- ggplot(M09_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M09_sgnf)

#save the plot
# genus: 3
ggsave("M09_sgnf.png", plot_M09_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M09_sgnf.svg", plot_M09_sgnf, width = 12.5, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M10
#plot
plot_M10_sgnf <- ggplot(M10_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M10_sgnf)

#save the plot
# genus: 2
ggsave("M10_sgnf.png", plot_M10_sgnf, width = 10, height = 15, units = "cm")
ggsave("M10_sgnf.svg", plot_M10_sgnf, width = 10, height = 15, units = "cm")


### Plot, print and save using ggplot2
## M11
#plot
plot_M11_sgnf <- ggplot(M11_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M11_sgnf)

#save the plot
# genus: 3
ggsave("M11_sgnf.png", plot_M11_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M11_sgnf.svg", plot_M11_sgnf, width = 12.5, height = 15, units = "cm")


## M12
#plot
plot_M12_sgnf <- ggplot(M12_sgnf, aes(x = genus, y = value, color = feeding)) +
  geom_boxplot(width = 0.5) +
  labs(title = "Boxplot Aldex2 diff abundant genera between feeding types", x = "genus", y = "clr_abundance") +
  geom_point(position = position_dodge(width = 0.5), size = 1)+
  theme_minimal()

#print
print(plot_M12_sgnf)

#save the plot
# genus: 3
ggsave("M12_sgnf.png", plot_M12_sgnf, width = 12.5, height = 15, units = "cm")
ggsave("M12_sgnf.svg", plot_M12_sgnf, width = 12.5, height = 15, units = "cm")


