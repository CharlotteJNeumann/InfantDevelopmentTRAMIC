#################################################################################

################################# Boxplots ######################################

#################################################################################




########################################################################
###################### ORAL ############################################
########################################################################

# Set working directory

setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/13_beta_diversity_and_Permanova/03_longitudinal/03_vizualization/")

# Load data
oral <- read.csv("LME_raw-data.csv", header = TRUE, sep = ",", row.names = 1)

# Exclude NAs in breastfeeding
oral2 <- subset(oral, breastfeeding == "x" | breastfeeding == "n")

# Ensure timepoint is a factor
oral2$timepoint <- as.factor(oral2$timepoint)

# Load ggplot2 library
library(ggplot2)

# Create the plot
distance <- ggplot(oral2, aes(x = timepoint, y = Distance, fill = breastfeeding)) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(width = 0.8), size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for better readability
print(distance)


########################################################################
###################### STOOL ###########################################
########################################################################

# Set working directory

setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/02_infants_development_stool_bacteria/13_beta_diversity_and_Permanova/03_longitudinal_LMR/03_vizualization/")

# Load data
stool <- read.csv("LME_raw-data.csv", header = TRUE, sep = ",", row.names = 1)

# Exclude NAs in breastfeeding
stool2 <- subset(stool, breastfeeding == "x" | breastfeeding == "n")

# Ensure timepoint is a factor
stool2$timepoint <- as.factor(stool2$timepoint)

# Load ggplot2 library
library(ggplot2)

# Create the plot
distance <- ggplot(stool2, aes(x = timepoint, y = Distance, fill = breastfeeding)) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(width = 0.8), size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for better readability
print(distance)
