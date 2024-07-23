###### stress_centrality_line_chart ########

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/01_infants_development_oral_bacteria/09_network/04_R_line_graph_centrality")

# Load required libraries
library(ggplot2)
library(tidyr)

#########################################
############# oral ######################
#########################################

############## BF ########################

# Read data from CSV
data <- read.csv("02_input/oral_BF_stress_centrality.csv", stringsAsFactors = FALSE)

# Convert data to long format
data_long <- gather(data, Bacteria, Value, -tp)


# Plot with dots
ggplot(data_long, aes(x = tp, y = Value, color = Bacteria, group = Bacteria)) +
  geom_line() +
  geom_point() +  # Add points
  labs(title = "oral BF stress centrality over time",
       y = "stress centrality") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))


############## NBF ########################

# Read data from CSV
data <- read.csv("02_input/oral_NBF_stress_centrality.csv", stringsAsFactors = FALSE)

# Convert data to long format
data_long <- gather(data, Bacteria, Value, -tp)


# Plot with dots
ggplot(data_long, aes(x = tp, y = Value, color = Bacteria, group = Bacteria)) +
  geom_line() +
  geom_point() +  # Add points
  labs(title = "oral NBF stress centrality over time",
       y = "stress centrality") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))




#########################################
############# stool ######################
#########################################

# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/02_infants_development_stool_bacteria/09_Network/04_R_line_graph_centrality")


############## BF ########################

# Read data from CSV
data <- read.csv("02_input/stool_BF_stress_centrality.csv", stringsAsFactors = FALSE)

# Convert data to long format
data_long <- gather(data, Bacteria, Value, -tp)


# Plot with dots
ggplot(data_long, aes(x = tp, y = Value, color = Bacteria, group = Bacteria)) +
  geom_line() +
  geom_point() +  # Add points
  labs(title = "stool BF stress centrality over time",
       y = "stress centrality") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))


############## NBF ########################

# Read data from CSV
data <- read.csv("02_input/stool_NBF_stress_centrality.csv", stringsAsFactors = FALSE)

# Convert data to long format
data_long <- gather(data, Bacteria, Value, -tp)


# Plot with dots
ggplot(data_long, aes(x = tp, y = Value, color = Bacteria, group = Bacteria)) +
  geom_line() +
  geom_point() +  # Add points
  labs(title = "stool NBF stress centrality over time",
       y = "stress centrality") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
