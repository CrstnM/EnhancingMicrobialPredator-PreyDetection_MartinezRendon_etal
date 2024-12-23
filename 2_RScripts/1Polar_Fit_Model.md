# Joint Species Distribution Models in R

Cristina Martínez Rendón  
04-12-2024  


Joint Species Distribution Modelling (JSDM) has emerged as a powerful and increasingly popular statistical approach for analyzing complex data in community ecology. Among the tools available, Hierarchical Modeling of Species Communities (HMSC) stands out as a versatile framework that applies the principles of generalized linear models (GLMs) within a Bayesian inference approach. HMSC allows the integration of community ecology data with data on environmental covariates, species traits, phylogenetic relationships, and the spatio-temporal context of the study, providing predictive insights into community assembly processes from non-manipulative observational data of species communities (Tikhonov et al. 2020) 
 
**R version:** 4.3.0 (21-04-2023)  

**Packages**  

``` r
# Install packages
# if (!require("Hmsc")) install.packages("Hmsc")
# if (!require("BayesLogit")) install.packages("BayesLogit")
# if (!require("corrplot")) install.packages("corrplot")
# if (!require("tidyverse")) install.packages("tidyverse")
# if (!require("RColorBrewer")) install.packages("RColorBrewer")
# if (!require("viridis")) install.packages("viridis")
# if (!require("plyr")) install.packages("plyr")
# if (!require("abind")) install.packages("abind")

# Load packages
library(devtools)
library(BayesLogit)
library(Hmsc)
library(corrplot)
library(tidyverse)
library(RColorBrewer)
library(coda)
library(viridis)
library(plyr)
library(abind)

rm(list = ls())

# Set directories
setwd("~/R_Projects/ArcticAntarctica/HMSC")

localDir = "."
dataDir = file.path(localDir, "Input_data")
     #dataDir = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/Input_data')
ModelDir = file.path(localDir, "models")
    #ModelDir = file.path('/home/alle/WORKING_DIR/Cristina/HMSC/models')
``` 

## 0. Defining the models
Sequencing data sets are characterized by a high frequency of zeros, indicating many species absences, which necessitated the use of a **hurdle model**. This modeling approach has two components: one for presence-absence and another for abundance conditioned on presence. (Following Odriozola et al. (2021)).

- Presence-Absence Modeling: The data were converted into a binary format, preserving zeros as they were and converting all nonzero values to one. This binary data set was analyzed using a binomial model with a probit link function, applied separately to each column (representing individual OTUs).

- Abundance Conditional on Presence: For the second component, all zeros in the data were replaced with missing values, while nonzero values retained their original measurements. This modified data set, representing scaled abundances (standardized to zero mean and unit variance), was modeled using the Possion lognormal distribution. For more information visit Ovaskainen & Abrego (2020). 


## 1. Read species data sets and wrangle them
``` r
  # Cercozoa data
  Cercozoa <- read.csv(file.path(dataDir,"CountData_Cercozoa.csv"), header = TRUE, row.names = 1, sep=";")
  colnames(Cercozoa) <- gsub("^X", "", colnames(Cercozoa))
  Cerco <- as.data.frame(t(Cercozoa))
  nyCer <- dim(Cerco)[1]
  
  # Green Algae data
  GreenAlgae <- read.csv(file.path(dataDir,"CountData_GreenAlgae.csv"), header = TRUE, row.names = 1, sep=";")
  colnames(GreenAlgae) <- gsub("^X", "", colnames(GreenAlgae))
  GAlgae <- as.data.frame(t(GreenAlgae))
  nyGAlgae <- dim(GAlgae)[1]
  
  # Ochrophyte data
  Ochrophytes <- read.csv(file.path(dataDir,"CountData_Diatoms.csv"), header = TRUE, row.names = 1, sep=";")
  colnames(Ochrophytes) <- gsub("^X", "", colnames(Ochrophytes))
  Ochro <- as.data.frame(t(Ochrophytes))
  nyOchro <- dim(Ochro)[1]
``` 
    

## 2. Prevalence filtering  
I included in the analyses only those OTUs with a prevalence of at least 10% among sampling units, consistent with the approach used for network analyses in FlashWeave. This threshold was chosen because data on rare species lack sufficient informativeness to support the fitting of reliable species-specific models (Ovaskainen & Abrego, 2020).
  
``` r  
  threshold.prev = 0.1
  
  #Cercozoa
  Cerco.pa<-ifelse(Cerco>0,1,0)
  Cerco.rel<-Cerco/rep(rowSums(Cerco),times=length(Cerco))
  cond1=!(colSums(Cerco.pa)<= threshold.prev*nyCer)
  Cercodata <- Cerco[,cond1]
  dim(Cercodata)

  #Green Algae
  GAlgae.pa<-ifelse(GAlgae>0,1,0)
  GAlgae.rel<-GAlgae/rep(rowSums(GAlgae),times=length(GAlgae))
  cond1=!(colSums(GAlgae.pa)<= threshold.prev*nyGAlgae)
  GAlgaedata <- GAlgae[,cond1]
  dim(GAlgaedata) 
  
  #Ochrophytes
  Ochro.pa<-ifelse(Ochro>0,1,0)
  Ochro.rel<-Ochro/rep(rowSums(Ochro),times=length(Ochro))
  cond1=!(colSums(Ochro.pa)<= threshold.prev*nyOchro)
  Ochrodata <- Ochro[,cond1]
  dim(Ochrodata)

  # Merge the three data sets
  Merged <- as.data.frame(cbind(Cercodata, GAlgaedata, Ochrodata))
  Y <- as.matrix(cbind(Cercodata, GAlgaedata, Ochrodata))
  rownames(Y) <- NULL
  Y <- as.matrix(Y)
  ny = dim(Y)[1]
  ns = dim(Y)[2]
```  

## 3. Read and log-transform environmental data (if necessary)
``` r
  Env <- read.csv(file.path(dataDir,"Polar_env.csv"), header = TRUE, sep=";")
  Env$sample_code <- as.factor(Env$sample_code)
  Env$set <- as.factor(Env$set)  
  Env$site <- as.factor(Env$site)  
  
  plot(Env)
  hist(Env$mean_N_per100)
  hist(Env$mean_C_per100)
  hist(Env$P_gperkg)
  Env$mean_N_per100 <- log(Env$mean_N_per100)
  Env$mean_C_per100 <- log(Env$mean_C_per100)
  
  # Depths:
  CercoDepth<-log(rowSums(Cerco))
  GAlgaeDepth<-log(rowSums(GAlgae))
  OchroDepth<-log(rowSums(Ochro))
  
  hist(CercoDepth)
  hist(GAlgaeDepth)
  hist(OchroDepth)
  
  XData <- data.frame(Env[,c(1,2,3,6,7,8,10)],CercoDepth,GAlgaeDepth,OchroDepth)
  
```   
## 4. Add trait data (TrData matrix)
``` r  
  Cerco_traits <- read.csv(file.path(dataDir,"Taxonomy_Cercozoa.csv"), header = TRUE, sep=";")
  GreenA_traits <- read.csv(file.path(dataDir,"Taxonomy_GreenAlgae.csv"), header = TRUE, sep=";")
  Ochro_traits <- read.csv(file.path(dataDir,"Taxonomy_Diatoms.csv"), header = TRUE, sep=";")
    all_traits <- rbind(Cerco_traits, GreenA_traits, Ochro_traits)
  
  Merged_t <- as.data.frame(t(Merged))
  Merged_t$OTU_Sp <- rownames(Merged_t)
  Merged_t <- Merged_t[, "OTU_Sp", drop = FALSE]
  Merged_t$Original_Order <- seq_len(nrow(Merged_t))
``` 
Match and merge the nutrition and phylogeny data into Merged_t
``` r
  Merged_t <- merge(Merged_t, 
                    all_traits[, c("OTU_Sp", "nutrition", "Kingdom", "Supergroup", 
                                   "Phylum", "Class", "Order", "Family", "Genus")],
                    by = "OTU_Sp", 
                    all.x = TRUE)
  Merged_t <- Merged_t %>% arrange(Original_Order)
  Merged_t$Original_Order <- NULL
```   
I saved the taxa order for use it to generate the species networks figures
``` r
  desired_order <- Merged_t$OTU_Sp
  write.csv(desired_order, file = file.path(dataDir, "desired_order.csv"), row.names = TRUE)
  
  TrData <- Merged_t
  rownames(TrData) <- TrData$OTU_Sp
  TrData$OTU_Sp <- NULL
  TrData <- TrData %>% 
    mutate(across(where(is.character), as.factor))
  TrData <- droplevels(TrData) # Ensure that all factor levels in TrData are valid (i.e., no unused or "empty" levels)
  dim(TrData)
```   
I saved as well trait data to filter out bacterivores in a second version of the species networks
``` r
  write.csv(TrData, file = file.path(dataDir, "TrData_subset.csv"), row.names = TRUE)
``` 

## 5. Hmsc-specific data wrangling
``` r
  
  #Abundance matrix:
    # Handle zeros in the abundance matrix Y
  Yabu = Y  # Copy of the merged data
  Yabu[Y == 0] = NA  # Replace zeros with NA
  Yabu = log(Yabu)  # Log-transform abundance
  Yabu[is.na(Yabu)] <- 0
  
  # Normalize / Standardize abundances for each species
  for (i in 1:ns) {
    Yabu[, i] = Yabu[, i] - mean(Yabu[, i], na.rm = TRUE)  # Center each species
    Yabu[, i] = Yabu[, i] / sd(Yabu[, i], na.rm = TRUE)    # Normalize by SD
  }
  
  # Quality check
  summary(Yabu)
  any(is.infinite(Yabu))  # Should return FALSE
  apply(Yabu, 2, mean, na.rm = TRUE)  # Should be near 0
  apply(Yabu, 2, sd, na.rm = TRUE)  # Should be near 1
  
  
  # Binary presence-absence matrix:
  Ypa = 1 * (Y > 0)
  
``` 
### Depth integration:
- Adjust depth for each taxonomic group using your TrData (merged taxonomy and nutrition data).
- Depth information should be dynamically selected or calculated for XData based on Phylum (e.g., Cercozoa, Green Algae, Ochrophytes).
``` r  
  Domain <- rep(rep(c("Cercozoa", "GreenAlgae", "Ochrophytes"), c(dim(Cercodata)[2], dim(GAlgaedata)[2], dim(Ochrodata)[2])))
  
  # Create a list for environment data adjusted for each species
  XDataList = list()
  
  # Loop over each species to adjust Depth based on their group (Domain)
    
      for (i in 1:length(Domain)) { # Inner loop for each entry in Domain
        tmp = XData
        
        # Dynamically assign depth based on Domain
        if (Domain[i] == "Cercozoa") {
          tmp$Depth = tmp$CercoDepth
        } else if (Domain[i] == "GreenAlgae") {
          tmp$Depth = tmp$GAlgaeDepth
        } else if (Domain[i] == "Ochrophytes") {
          tmp$Depth = tmp$OchroDepth
        }
        
        # Store the modified XData in XDataList
        XDataList[[i]] = tmp
      }
  
    
    
```  
## 6. MODEL SETUP!
  
### Study design (RANDOM EFFECTS)
- Random effects model variation at each levels that won't be taken as fixed effects, but rather as sources of variability. This is commonly done in hierarchical or mixed-effects models.
- The rL.x objects specify random effects. The setPriors() function determines the number of latent factors (nfMin and nfMax) that can explain random variation. This is where spatial or sample-level variability is captured, critical for community ecology models.
``` r
  studyDesign = data.frame(sample_code=XData$sample_code, site=XData$site)
    studyDesign$site = factor(studyDesign$site, levels = unique(studyDesign$site))
    studyDesign$sample_code = factor(studyDesign$sample_code, levels = unique(studyDesign$sample_code))
  
  # Set random effects
  rL.site = HmscRandomLevel(units = levels(studyDesign$site))
  rL.sample_code = HmscRandomLevel(units = levels(studyDesign$sample_code))
  
  
  # Formula(s) (FIXED EFFECTS)
  # Regression models for environmental covariates (predictor variables of species abundance, fixed effects)
  XFormula1 = ~ Depth # Null model has only the sequencing depth as the explanatory variable
  XFormula2 = ~ set + ph + mean_N_per100 + mean_C_per100 + P_gperkg + Depth # Full model
  
  # Regression model for traits
  TrFormula = ~ nutrition + Phylum
  
```  
#### MODEL = 1: NULL MODELS (XFormula1); MODEL = 2: ENVIRONMENTAL AND SPATIAL PREDICTORS (XFormula); 
  
CONSTRUCT THE MODELS 
``` r 
  create_hmsc_model <- function(Y, XFormula, distribution) {
    Hmsc(
      Y = Y,
      XData = XDataList, #XData
      XFormula = XFormula, 
      TrData = TrData,
      TrFormula = TrFormula,
      distr = distribution,
      studyDesign = studyDesign,
      ranLevels = list(site = rL.site, sample_code = rL.sample_code)
    )
  }
```   
Presence-absence models (Null and full). 
``` r
  modpa_null <- create_hmsc_model(Y = Ypa, XFormula = XFormula1, distribution = "probit")
  modpa_full <- create_hmsc_model(Y = Ypa, XFormula = XFormula2, distribution = "probit")
```   
  
Abundance models (Null and full)
``` r
  modabu_null <- create_hmsc_model(Y = Yabu, XFormula = XFormula1, distribution = "normal")
  modabu_full <- create_hmsc_model(Y = Yabu, XFormula = XFormula2, distribution = "normal")
```   
  
COMBINING AND SAVING MODELS
``` r  
  models <- list(
    modpa_null = modpa_null,
    modpa_full = modpa_full,
    modabu_null = modabu_null,
    modabu_full = modabu_full)

  names(models) = c(modpa_null, modpa_full, modabu_null, modabu_full)
  save(models, file = file.path(ModelDir, "unfitted_models.RData"))
```   
  
TESTING THAT MODELS FIT WITHOUT ERRORS. 
``` r  
  for(i in 1:length(models)){
    print(i)
    sampleMcmc(models[[i]],samples=2)
  }
  
```   
Sampling the model (Test) 
``` r  
  set.seed(850511)
  #dir.create("~/R_Projects/ArcticAntarctica/HMSC/Test")
  
  nParallel = NULL #Default: nParallel = nChains
  
  
  ## Run HMSC tests for all models 180 iterations. (For loop a nightmare, repetitive code follows.)
  samples = 20 # Posterior samples per chain after burn-in
  thin = 10    # Thinning, i.e., between how many samples one is kept
  transient = ceiling(samples*thin*0.4) # Discarded iterations, here are 80 iterations discarded
  nChains = 4
 

  t0 <- Sys.time()
  hmsc_modabu_null2 <- "~/R_Projects/ArcticAntarctica/HMSC/Test/hmsc_modabu_null2.Rda"
  if(!file.exists(hmsc_modabu_null2)){
    modabu_null = sampleMcmc(modabu_null, samples = samples, thin = thin,
                   transient = transient, #adaptNf = adaptNf,
                   nChains = nChains, nParallel = nChains, verbose = 2)
    print(Sys.time()-t0)
    save(modabu_null, file=hmsc_modabu_null2)
  }else{load(hmsc_modabu_null2)}

        t0 <- Sys.time()
        hmsc_modabu_full1 <- "~/R_Projects/ArcticAntarctica/HMSC/Test/hmsc_modabu_full1.Rda"
        if(!file.exists(hmsc_modabu_full1)){
          modabu_full = sampleMcmc(modabu_full, samples = samples, thin = thin,
                                   transient = transient, #adaptNf = adaptNf,
                                   nChains = nChains, nParallel = nChains, verbose = 2)
          print(Sys.time()-t0)
          save(modabu_full, file=hmsc_modabu_full1)
        }else{load(hmsc_modabu_full1)}

# Run HMSC for presence-absence data with 140 iterations.

  t0 <- Sys.time()
  hmsc_modpa_null1 <- "~/R_Projects/ArcticAntarctica/HMSC/Test/hmsc_modpa_null1.Rda"
  if(!file.exists(hmsc_modpa_null1)){
    modpa_null = sampleMcmc(modpa_null, samples = samples, thin = thin,
                             transient = transient, #adaptNf = adaptNf,
                             nChains = nChains, nParallel = nChains, verbose = 2)
    print(Sys.time()-t0)
    save(modpa_null, file=hmsc_modpa_null1)
  }else{load(hmsc_modpa_null1)}

        t0 <- Sys.time()
        hmsc_modpa_full1 <- "~/R_Projects/ArcticAntarctica/HMSC/Test/hmsc_modpa_full1.Rda"
        if(!file.exists(hmsc_modpa_full1)){
          modpa_full = sampleMcmc(modpa_full, samples = samples, thin = thin,
                                   transient = transient, #adaptNf = adaptNf,
                                   nChains = nChains, nParallel = nChains, verbose = 2)
          print(Sys.time()-t0)
          save(modpa_full, file=hmsc_modpa_full1)
        }else{load(hmsc_modpa_full1)}


    gc() # (To reduce the chance of memory bloat after many iterations, especially with large Bayesian models).
```  
    

## 7. Scale up the sampling process
Run HMSC test for the four models

``` r     
  # Sampling parameters
    thin = 100          # Keep every 100th iteration a posterior sample.      
    samples = 1000      # This is the target number of posterior samples *per chain* after thinning.  
    transient = ceiling(samples*thin*0.5) #It does 50k iterations and discard them (burn-in) before taking samples.
    nChains = 4         # Run 4 MCMC chains in parallel
    #ModelDir <- "~/R_Projects/ArcticAntarctica/HMSC/Results"
        #ModelDir <-'/home/alle/WORKING_DIR/Cristina/HMSC/Results'
    
    # Create output directory if not exists
    if (!dir.exists(ModelDir)) dir.create(ModelDir, recursive = TRUE)
    
    set.seed(850511)
    
    # Loop through models to run sampling and save results
    for (model_name in names(models)) {
      cat("Running model:", model_name, "\n")
      
      # Start the timer
      ptm = proc.time()
      
      # Sample MCMC
      sampled_model <- sampleMcmc(
        models[[model_name]],
        samples = samples,
        thin = thin,
        transient = transient,
        nChains = nChains,
        nParallel = nChains,
        verbose = 2
      )
      
      # Calculate computational time
      computational.time =  proc.time() - ptm
      

      # Generate filenames
      model_filename <- file.path(ModelDir, paste0(model_name, "_thin_", thin, "_samples_", samples, "_chains_", nChains, ".rds"))
      time_filename <- file.path(ModelDir, paste0(model_name, "_thin_", thin, "_samples_", samples, "_chains_", nChains, "_comptime.txt"))
      
      # Save sampled model and computational time
      saveRDS(sampled_model, file = model_filename)
      write.table(as.matrix(computational.time), file = time_filename, row.names = FALSE, col.names = FALSE)
    }
    
    cat("All models processed and saved.\n")

    gc()  
``` 
 


References:
- Tikhonov G, Opedal ØH, Abrego N, Lehikoinen A, de Jonge MMJ, Oksanen J, et al. Joint species distribution modelling with the r-package Hmsc. Methods Ecol Evol. 2020;11:442–7. 
- Odriozola I, Abrego N, TlÁskal V, ZrůstovÁ P, Morais D, Větrovský T, et al. Fungal Communities Are Important Determinants of Bacterial Community Composition in Deadwood. Shade A, editor. mSystems. 2021;6:e01017-20. 
- Ovaskainen O, Abrego N. Joint Species Distribution Modelling: With Applications in R. Cambridge: Cambridge University Press; 2020. Available from: https://www.cambridge.org/core/product/0D9FA93EA1DD408332A17266449668B3