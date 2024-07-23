# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/03_infants_development_oral_archaea/11_statistics")

# read file
oral_archaea <- read.csv("inf_dev_oral_archaea_pres_abs.csv", header = TRUE)


####################################################
##################### oral #########################
################### feeding ########################
####################################################

##### Methanobrevibacter

# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobrevibacter, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 0.53947368 1.00000000 1.00000000 0.08235294 1.00000000 0.48928571 1.00000000 0.62669683 1.00000000
#[10] 1.00000000 1.00000000 1.00000000 0.39560440         NA


### unclassified_Woesearchaeales

# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#  [1] 0.5597910 1.0000000 0.6328671 0.4264706 0.4725275 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
# [11] 0.5300699 0.4909091 0.2747253        NA


### Candidatus_Nitrosotenuis
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Candidatus_Nitrosotenuis, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]  1  1 NA NA NA  1  1 NA NA NA NA NA NA NA


### Methanobacterium
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobacterium, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#  [1] 1.0000000 0.2222222        NA 1.0000000        NA        NA        NA        NA        NA 1.0000000
# [11] 1.0000000        NA 1.0000000        NA


### Methanosphaera
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanosphaera, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#  [1] NA NA NA  1 NA NA NA  1 NA NA NA NA NA NA



### unclassified_Aenigmarchaeales
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Aenigmarchaeales, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] NA NA NA NA NA NA NA NA NA NA NA  1 NA NA


### unclassified_Nitrososphaeraceae
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Nitrososphaeraceae, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]     NA 1.0000     NA     NA     NA 0.1875     NA     NA 1.0000     NA 1.0000 1.0000 1.0000     NA



### unclassified_Woesearchaeales_GW2011_GWC1_47_15
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales_GW2011_GWC1_47_15, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#  [1]        NA        NA        NA 0.2352941 1.0000000        NA        NA        NA        NA        NA
# [11]        NA        NA        NA        NA



### unclassified_Woesearchaeales_SCGC_AAA011-D5
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales_SCGC_AAA011_D5, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] NA NA NA  1 NA  1 NA NA NA NA NA NA NA NA



#####################################################
################## oral #############################
################# birthmode #########################
#####################################################



##### Methanobrevibacter

# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobrevibacter, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 1.0000000 1.0000000 0.5879121 1.0000000 1.0000000 0.5764706 0.1515152 0.6266968 1.0000000 1.0000000
#[11] 1.0000000 0.1818182 1.0000000 1.0000000


### unclassified_Woesearchaeales

# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 0.6026832 1.0000000 0.6328671 0.4705882 0.2773893 0.4705882 1.0000000 1.0000000 0.5151515 0.4736842
#[11] 1.0000000 0.5454545 1.0000000 0.303405


### Candidatus_Nitrosotenuis
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Candidatus_Nitrosotenuis, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 1.0000000 0.4771242        NA        NA        NA 0.4705882 0.4166667        NA        NA        NA
#[11]        NA        NA        NA        NA


### Methanobacterium
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobacterium, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 0.4000000 0.4444444        NA 0.2941176        NA        NA        NA        NA        NA 1.0000000
#[11] 1.0000000        NA 1.0000000 0.3698500


### Methanosphaera
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanosphaera, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1]        NA        NA        NA 1.0000000        NA        NA        NA 1.0000000        NA        NA
#[11]        NA        NA        NA 0.4736842



### unclassified_Aenigmarchaeales
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Aenigmarchaeales, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
#[11]        NA 0.4545455        NA 0.4736842


### unclassified_Nitrososphaeraceae
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Nitrososphaeraceae, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1]        NA 0.4444444        NA        NA        NA 1.0000000        NA        NA 1.0000000        NA
#[11] 1.0000000 0.4545455 0.3571429        NA


### unclassified_Woesearchaeales_GW2011_GWC1_47_15
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales_GW2011_GWC1_47_15, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] NA NA NA  1  1 NA NA NA NA NA NA NA NA  1



### unclassified_Woesearchaeales_SCGC_AAA011-D5
# Get unique time points
unique_time_points <- unique(oral_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- oral_archaea[oral_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales_SCGC_AAA011_D5, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1]        NA        NA        NA 1.0000000        NA 0.4705882        NA        NA        NA        NA
#[11]        NA        NA        NA 0.2105263


###################################################
################# STOOL ###########################


# set working directory
setwd("C:/Users/o_neumannj/Nextcloud/Charlotte/TRAMIC/03_infant_development/Amplikon_Analysis/04_infants_development_stool_archaea/12_statistics")

# read file
stool_archaea <- read.csv("inf_dev_stool_archaea_pres_abs.csv", header = TRUE)


####################################################
##################### stool #########################
################### feeding ########################
####################################################

##### Methanobrevibacter

# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobrevibacter, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 0.40000000         NA 0.14285714         NA 0.52380952 0.14285714 1.00000000         NA 0.06060606
#[10]         NA 0.14285714 1.00000000         NA 1.00000000 1.00000000


### Halococcus

# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Halococcus, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1]        NA        NA        NA 1.0000000 1.0000000 1.0000000        NA        NA 0.1818182        NA
#[11] 1.0000000        NA        NA 1.0000000        NA


### Methanobacterium
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobacterium, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] NA NA NA NA NA NA  1 NA  1 NA NA NA NA NA NA


### Methanosarcina
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanosarcina, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 0.3333333        NA        NA        NA 1.0000000        NA        NA        NA 0.4545455        NA
# [11]        NA        NA        NA 1.0000000 1.0000000


### Methanosphaera
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanosphaera, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 0.3333333 1.0000000 1.0000000        NA 0.4444444 1.0000000 0.1538462        NA 0.4545455        NA
# [11] 0.1071429        NA        NA 0.1428571 1.0000000 



### Thermococcus
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Thermococcus, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]        NA        NA        NA        NA        NA        NA        NA        NA 0.4545455        NA
#[11] 1.0000000        NA        NA        NA        NA


### unclassified_Nitrososphaeraceae
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Nitrososphaeraceae, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 0.46666667         NA 0.42857143 1.00000000 0.16666667 0.48571429 1.00000000         NA 0.06060606
# [10] 1.00000000 1.00000000 1.00000000         NA         NA 1.00000000


### unclassified_Woesearchaeales
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1]        NA        NA        NA        NA        NA        NA 1.0000000        NA 0.4545455        NA
#[11]        NA        NA        NA        NA        NA


### Candidatus_Nitrocosmicus
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Candidatus_Nitrocosmicus, subset_data$breastfeeding)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 1.0 0.4  NA  NA 1.0  NA  NA  NA 1.0 1.0 1.0 1.0  NA 1.0 0.4



#####################################################
################## stool #############################
################# birthmode #########################
#####################################################



##### Methanobrevibacter

# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobrevibacter, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 1.0000000        NA 1.0000000        NA 0.1939394 1.0000000 0.2861305 1.0000000 1.0000000        NA
# [11] 1.0000000 0.2857143 1.0000000 1.0000000 0.4000000


### Candidatus_Nitrocosmicus

# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Candidatus_Nitrocosmicus, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
#[1] 1.0000000 0.4000000        NA        NA 0.3636364        NA        NA        NA 0.4545455 1.0000000
# [11] 1.0000000 0.2857143        NA 0.2857143 1.0000000


### Halococcus
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Halococcus, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]        NA        NA        NA 1.0000000 1.0000000 0.4285714        NA        NA 0.4545455        NA
# [11] 0.3333333        NA        NA 1.0000000        NA


### Methanobacterium
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanobacterium, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]        NA        NA        NA        NA        NA        NA 1.0000000 0.2305311 0.4545455        NA
# [11]        NA        NA 0.5920746        NA        NA


### Methanosphaera
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanosphaera, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 1.0000000 1.0000000 1.0000000        NA 1.0000000 1.0000000 0.4615385 0.4545455 0.1818182        NA
# [11] 0.5000000        NA 0.4615385 1.0000000 1.0000000



### Methanosarcina
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Methanosarcina, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]  1 NA NA NA  1 NA NA  1  1 NA NA NA NA  1  1


### Thermococcus
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$Thermococcus, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]        NA        NA        NA        NA        NA        NA        NA        NA 1.0000000        NA
# [11] 0.3333333        NA        NA        NA        NA


### unclassified_Nitrososphaeraceae
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Nitrososphaeraceae, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1] 1.0000000        NA 1.0000000 1.0000000 1.0000000 0.4857143 1.0000000        NA 1.0000000 1.0000000
# [11] 1.0000000 0.2857143 0.4615385        NA 0.4000000



### unclassified_Woesearchaeales
# Get unique time points
unique_time_points <- unique(stool_archaea$timepoint)

# Initialize a vector to store p-values
p_values <- numeric(length(unique_time_points))

# Loop through each time point
for (i in seq_along(unique_time_points)) {
  # Subset the data for the current time point
  subset_data <- stool_archaea[stool_archaea$timepoint == unique_time_points[i], ]
  
  # Create a contingency table
  contingency_table <- table(subset_data$unclassified_Woesearchaeales, subset_data$birthmode)
  
  # Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Store the p-value for the current time point
    p_values[i] <- fisher_test$p.value
  } else {
    # If the contingency table doesn't meet the condition, set p-value to NA
    p_values[i] <- NA
  }
}

# Print the p-values for each time point
print(p_values)
# [1]        NA        NA        NA        NA        NA        NA 0.4615385 1.0000000 1.0000000        NA
# [11]        NA        NA        NA        NA        NA