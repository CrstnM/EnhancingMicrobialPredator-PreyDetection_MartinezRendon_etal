
### Model diagnostics: To test the model's convergence 
# Post-processing of model results

# Cristina Martínez Rendón
# 05-12-2024

# Load packages
library(Hmsc)
library(ggtext)
library(knitr)
library(plyr)
library(tibble)
library(ggplot2)
library(coda)


# Directories:
setwd("~/R_Projects/ArcticAntarctica/HMSC")

localDir <- "."

ModelDir <- file.path(localDir, "models")
#ModelDir <- file.path('/home/alle/WORKING_DIR/Cristina/HMSC/models')

ServerDir <- file.path(localDir, "Server_results")
#ServerDir <- file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Results')

DiagDir <- file.path(localDir, "Diagnostics_results")
#DiagDir <- file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Diagnostics_results')
if (!dir.exists(DiagDir)) {
  dir.create(DiagDir, recursive = TRUE)
}

# Load data
modabu_null <- readRDS(file.path(ServerDir, "modabu_null_thin_100_samples_1000_chains_4.rds"))
modabu_full <- readRDS(file.path(ServerDir, "modabu_full_thin_100_samples_1000_chains_4.rds"))
modpa_null <- readRDS(file.path(ServerDir, "modpa_null_thin_100_samples_1000_chains_4.rds"))
modpa_full <- readRDS(file.path(ServerDir, "modpa_full_thin_100_samples_1000_chains_4.rds"))

# Models List
ranmodels_list <- list(modabu_null = modabu_null, 
                       modabu_full = modabu_full, 
                       modpa_null = modpa_null, 
                       modpa_full = modpa_full)


### 1. Diagnostics 

# - Effective Sample Size (ESS): Calculates the effective sample size of the Beta, Gamma, and V parameters to assess how well the MCMC chains are mixing.
# - Gelman Diagnostic (PSRF): Uses the Gelman-Rubin diagnostic to evaluate convergence, where values close to 1 indicate good convergence.

for (model_name in names(ranmodels_list)) {
  model <- ranmodels_list[[model_name]]
  mpost <- convertToCodaObject(model)
  
  #Examine MCMC and effective sample size
  es.beta = effectiveSize(mpost$Beta)
  ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  
  es.gamma = effectiveSize(mpost$Gamma)
  ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
  
  es.V = effectiveSize(mpost$V)
  ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
  
  # Effective Sample Size Diagnostics
  ess_beta_tibble <- tibble(ess_beta = as.numeric(es.beta))
  
  # Gelman diagnostic: should be <1.001
  psrf_beta_tibble <- tibble(psrf_beta = as.numeric(ge.beta[, "Point est."]))
  
  # Visualization
    # Effective Sample Size (ESS) plot
  ess_plot <- ggplot(ess_beta_tibble, aes(x = ess_beta)) +
    geom_histogram(binwidth = 10, color = "black", fill = "steelblue") +
    xlab("Effective Sample Size") +
    ggtitle(paste("ESS for", model_name)) +
    theme_minimal()
  
    # Gelman Diagnostic (PSRF) plot
  psrf_plot <- ggplot(psrf_beta_tibble, aes(x = psrf_beta)) +
    geom_histogram(binwidth = 0.01, color = "black", fill = "coral") +
    xlab("Gelman Diagnostic") +
    ggtitle(paste("PSRF for", model_name)) +
    theme_minimal()
  
  # Save plots
  ess_filename <- file.path(DiagDir, paste0(model_name, "_ess_plot.tiff"))
  psrf_filename <- file.path(DiagDir, paste0(model_name, "_psrf_plot.tiff"))
  ggsave(ess_plot, filename = ess_filename, width = 7, height = 5, bg = "white")
  ggsave(psrf_plot, filename = psrf_filename, width = 7, height = 5, bg = "white") 
  
   # Save all diagnostics as RDA
   save_filename <- file.path(DiagDir, paste0(model_name, "_diagnostics.rda"))
   
}
  


### 2.  Predictions and model fit
# The code evaluates the model fit and predictive performance of the HMSC models. It computes predicted values for each model using computePredictedValues and evaluates explanatory power (e.g., R² for abundance models, AUC for presence-absence models) and predictive power via cross-validation. 

# Initialize result storage
result_list <- list()

# Initialize an empty summary table
summary_table <- tibble(
  model = character(),
  type = character(),
  explanatory_power = numeric(),
  predictive_power = numeric()
)

# Loop through models for predictions and model fit
for (model_name in names(ranmodels_list)) {
  model <- ranmodels_list[[model_name]]
  
  # Compute predictions
  cat("Computing predictions for model:", model_name, "\n")
  preds <- computePredictedValues(model)
  
  # Evaluate model fit
  cat("Evaluating model fit for model:", model_name, "\n")
  MF <- evaluateModelFit(hM = model, predY = preds)
  
  # Save predictions and model fit results
  result_list[[model_name]] <- list(predictions = preds, model_fit = MF)
  
  # Save results as .rds
  saveRDS(list(preds = preds, MF = MF), file = file.path(DiagDir, paste0(model_name, "_predictions_and_fit.rds")))
  
  # Determine model type (abundance or presence-absence)
  type <- ifelse(grepl("abu", model_name), "abundance", "presence-absence")
  
  # Add row to the summary table
  summary_table <- summary_table %>% add_row(
    model = model_name,
    type = type,
    explanatory_power = if (type == "presence-absence") mean(MF$AUC, na.rm = TRUE) else mean(MF$R2, na.rm = TRUE),
    predictive_power = if (type == "presence-absence") mean(MF$AUC.CV, na.rm = TRUE) else mean(MF$R2.CV, na.rm = TRUE)
  )
}

# Save summary table
write.csv(summary_table, file = file.path(DiagDir, "summary_table.csv"), row.names = FALSE)


