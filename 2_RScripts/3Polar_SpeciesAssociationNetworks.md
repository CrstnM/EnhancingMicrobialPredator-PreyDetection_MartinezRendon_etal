# Post-processing of model results - Species association networks with HMSC models

Cristina Martínez Rendón  
06-12-2024  


**R version:** 4.3.0 (21-04-2023)  

**Packages**  

``` r
# Load packages
library(Hmsc)
library(ggtext)
library(knitr)
library(dplyr)
library(ggplot2)
library(ggcorrplot)

rm(list = ls())

# Set directories:
setwd("~/R_Projects/ArcticAntarctica/HMSC")

localDir = "."

dataDir = file.path(localDir, "Input_data")
   #dataDir = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Input_data')

ServerDir = file.path(localDir, "Server_results")
   #ServerDir = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Results')

NetworksDir = file.path(localDir, "Networks")
   #NetworksDir = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Networks')
if (!dir.exists(NetworksDir)) {
  dir.create(NetworksDir, recursive = TRUE)
}

Figures = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Figures')
if (!dir.exists(Figures)) {
  dir.create(Figures, recursive = TRUE)
}

toPlot_dir = file.path(localDir, "Network_matrices")
  #toPlot_dir = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Network_matrices')
if (!dir.exists(toPlot_dir)) {
  dir.create(toPlot_dir, recursive = TRUE)
}
```
Load data
``` r 
modabu_null <- readRDS(file.path(ServerDir, "modabu_null_thin_100_samples_1000_chains_4.rds"))
modabu_full <- readRDS(file.path(ServerDir, "modabu_full_thin_100_samples_1000_chains_4.rds"))
modpa_null <- readRDS(file.path(ServerDir, "modpa_null_thin_100_samples_1000_chains_4.rds"))
modpa_full <- readRDS(file.path(ServerDir, "modpa_full_thin_100_samples_1000_chains_4.rds"))

# Models List
ranmodels_list <- list(modabu_null = modabu_null, 
                       modabu_full = modabu_full, 
                       modpa_null = modpa_null, 
                       modpa_full = modpa_full)
```

## 1. Read taxa order to generate the species networks figures
```r 
desired_order <- read.csv(
  file = file.path(dataDir, "desired_order.csv"))$x
```

Loop through models for association computations. This loop computes the complete networks.
```r
for (model_name in names(ranmodels_list)) {
  model <- ranmodels_list[[model_name]]
  
  # Compute associations
  OmegaCor <- computeAssociations(model)
  supportLevel <- 0.89
  
  # Save matrices
  mean_filename <- file.path(NetworksDir, paste0("OmegaCor_mean_", model_name, ".csv"))
  support_filename <- file.path(NetworksDir, paste0("OmegaCor_support_", model_name, ".csv"))
  write.csv(OmegaCor[[1]]$mean, file = mean_filename)
  write.csv(OmegaCor[[1]]$support, file = support_filename)
#}
 
  # Generate matrices for mean and support
  hmdf_mean <- OmegaCor[[1]]$mean %>%
    as.matrix()
  hmdf_support <- OmegaCor[[1]]$support %>%
    as.matrix()
  
  # Reorder rows and columns to match the desired order
  hmdf_mean <- hmdf_mean[desired_order, desired_order]
  hmdf_support <- hmdf_support[desired_order, desired_order]
  
  # Filter associations based on supportLevel
  toPlot <- ((hmdf_support > supportLevel) + 
               (hmdf_support < (1 - supportLevel)) > 0) * hmdf_mean
  
  # Generate association plot
  omega_plot <- ggcorrplot::ggcorrplot(
    toPlot, 
    type = "lower", 
    hc.order = FALSE,  # Disable hierarchical clustering to preserve order
    title = paste("Species Associations for", model_name)
  ) +
    theme(
      axis.text.x = element_text(size = 4),    # X-axis text size
      axis.text.y = element_text(size = 4),    # Y-axis text size
      axis.title.x = element_text(size = 12), # X-axis title size
      axis.title.y = element_text(size = 12)  # Y-axis title size
    )
  
  # Save the plot
  plot_filename <- file.path("Figures", paste0("species_associations_", model_name, ".png"))
  ggsave(
    plot = omega_plot, 
    filename = plot_filename, 
    bg = "white", 
    width = 49, 
    height = 49
  )
}
```
## 2.  Trait-based filtering to exclude bacterivores, 
For a second version of the species networks, which includes only algivores and omnivores.
``` r 
Trait_subset <- read.csv(file = file.path(dataDir, "TrData_subset.csv"), header = TRUE, row.names = 1, sep=";") %>% 
  filter(nutrition %in% c("eukaryvore", "omnivore", "autotroph"))
  
desired_order <- rownames(Trait_subset) %>% 
  as.character(desired_order)

models_list <- list(
  modpa_null = computeAssociations(modpa_null),
  modpa_full = computeAssociations(modpa_full),
  modabu_null = computeAssociations(modabu_null),
  modabu_full = computeAssociations(modabu_full)
)

# Initialize a list to store plots and filtered matrices
plots <- list()

# Loop through each model
for (model_name in names(models_list)) {
  # Extract support and mean matrices
  OmegaCor <- models_list[[model_name]]
  hmdf_mean <- OmegaCor[[1]]$mean
  hmdf_support <- OmegaCor[[1]]$support
  
  # Subset rows and columns for both matrices
  hmdf_mean_subset <- hmdf_mean[desired_order, desired_order]
  hmdf_support_subset <- hmdf_support[desired_order, desired_order]
  
  # Filter associations based on support level
  supportLevel <- 0.89
  toPlot <- ((hmdf_support_subset > supportLevel) + 
               (hmdf_support_subset < (1 - supportLevel)) > 0) * hmdf_mean_subset
  
  # Save the toPlot matrix as a CSV file
  toPlot_filename <- file.path(toPlot_dir, paste0("toPlot_", model_name, ".csv"))
 ```
Here, I saved the matrix as a CSV file and later manually filtered it in Excel to retain only cross-phylum interactions, which represent the predator-prey interactions of interest.
``` 
  write.csv(toPlot, file = toPlot_filename)
  
#}
```
    
  # Generate association plot
``` r  
  omega_plot <- ggcorrplot::ggcorrplot(
    toPlot, 
    type = "lower", 
    hc.order = FALSE,  # Disable hierarchical clustering to preserve order
    title = paste("Species associations for", model_name,": Cercozoan algivores|omnivores with microalge")
  ) +
    theme(
      axis.text.x = element_text(size = 7),    # X-axis text size
      axis.text.y = element_text(size = 7),    # Y-axis text size
      axis.title.x = element_text(size = 16), # X-axis title size
      axis.title.y = element_text(size = 16),  # Y-axis title size
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)  # Title size, bold, and centered
    )
  
  # Save the plot
  plot_filename <- file.path(Figures, paste0("species_associations_", model_name, ".pdf"))
  ggsave(
    plot = omega_plot, 
    filename = plot_filename, 
    bg = "white", 
    width = 22, 
    height = 22
  )
  
  # Store plot in list for further use (optional)
  plots[[model_name]] <- omega_plot
}
  
``` 

## 3.  Trait-based filtering to exclude algivores and omnivores
A third version of the species networks.
``` r
Trait_subset2 <- read.csv(file = file.path(dataDir, "TrData_subset.csv"), header = TRUE, row.names = 1, sep=";") %>% 
  filter(nutrition %in% c("bacterivore", "autotroph"))

desired_order2 <- rownames(Trait_subset2)

# Initialize a list to store plots and filtered matrices
plots <- list()

# Loop through each model
for (model_name in names(models_list)) {
  # Extract support and mean matrices
  OmegaCor <- models_list[[model_name]]
  hmdf_mean <- OmegaCor[[1]]$mean
  hmdf_support <- OmegaCor[[1]]$support
  
  # Subset rows and columns for both matrices
  hmdf_mean_subset <- hmdf_mean[desired_order2, desired_order2]
  hmdf_support_subset <- hmdf_support[desired_order2, desired_order2]
  
  # Filter associations based on support level
  supportLevel <- 0.89
  toPlot <- ((hmdf_support_subset > supportLevel) + 
                              (hmdf_support_subset < (1 - supportLevel)) > 0) * hmdf_mean_subset
  
  # Save the toPlot matrix as a CSV file
  toPlot_filename <- file.path(toPlot_dir, paste0("toPlot_", model_name, ".csv"))
  write.csv(toPlot, file = toPlot_filename)
#}  
  # Generate association plot
  omega_plot <- ggcorrplot::ggcorrplot(
    toPlot, 
    type = "lower", 
    hc.order = FALSE,  # Disable hierarchical clustering to preserve order
    title = paste("Species associations for", model_name,": Cercozoan bacterivores with microalge")
  ) +
    theme(
      axis.text.x = element_text(size = 5),    # X-axis text size
      axis.text.y = element_text(size = 5),    # Y-axis text size
      axis.title.x = element_text(size = 16), # X-axis title size
      axis.title.y = element_text(size = 16),  # Y-axis title size
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)  # Title size, bold, and centered
    )
  
  # Save the plot
  plot_filename <- file.path(Figures, paste0("species_associations_", model_name, ".pdf"))
  ggsave(
    plot = omega_plot, 
    filename = plot_filename, 
    bg = "white", 
    width = 32, 
    height = 32
  )
  
  # Store plot in list for further use (optional)
  plots[[model_name]] <- omega_plot
}


### Count Zeros, Positive, and Negative Values

# Get a list of all CSV files in the directory
csv_files <- list.files(toPlot_dir, pattern = "*.csv", full.names = TRUE)

# Initialize a data frame to store the results
results <- data.frame(
  File = character(),
  Positive = integer(),
  Negative = integer(),
  Zeros = integer(),
  stringsAsFactors = FALSE
)

# Loop through each file and compute counts
for (file in csv_files) {
  # Read the CSV file
  matrix_data <- as.matrix(read.csv(file, row.names = 1))
  
  # Count positive, negative, and zero values
  positive_count <- sum(matrix_data > 0)
  negative_count <- sum(matrix_data < 0)
  zero_count <- sum(matrix_data == 0)
  
  # Store the results in the data frame
  results <- rbind(
    results,
    data.frame(
      File = basename(file),
      Positive = positive_count,
      Negative = negative_count,
      Zeros = zero_count,
      stringsAsFactors = FALSE
    )
  )
}

# Print the results
print(results)

# Optionally, save the results to a CSV file
write.csv(results, file = file.path(toPlot_dir, "toPlot_summary_bacterivores.csv"), row.names = FALSE)
```  
![Final species association networks after refining in Inkscape]
(../4_Figures/HMSC_Network.png)

