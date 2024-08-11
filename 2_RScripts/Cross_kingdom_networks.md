# Cross-kingdom network analysis  


Following the code of Jule Freudenthal  
Adapted by Cristina Martínez Rendón  
21-11-2023  

**R version:** 4.3.0 (21-04-2023)

**Packages**  
 
``` r
if (!require("compositions")) install.packages("compositions")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("vegan")) install.packages("vegan")

Load packages
library(compositions)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(vegan)
library(reshape2)


setwd("~/R_Projects/ArcticAntarctica/Polar_Networks")

rm(list = ls())
``` 

## 1. Data handling
I started by organizing the data for network analysis. Initially, I focused on the Arctic samples by filtering the metadata and loading the corresponding count and taxonomy tables for Cercozoa, Ochrophytes (here, Diatoms), and Green Algae. To reduce spurious putative interactions and based on NMDS analyses visualized from beta diversity dissimilarities, I conducted separate network analyses for each polar region: Arctic, Continental Antarctica, and the Antarctic Peninsula. This involved running the code independently for each dataset with minor adjustments, carefully marked with a ❄️ to ensure accuracy, to prevent errors during execution, and to avoid overwritting the generated files.
``` r  
# Create lists for count and taxonomy data
CountData <- list()
TaxonomyData <- list()

# Initialize empty list to store transformed data
Transposed_CountData <- list()

# Load metadata
metadata <- read.table("Metadata.csv", header = TRUE, sep = ";", dec = ".", row.names = 1)

# ❄️ Filter metadata
metadata <- metadata %>% 
  filter(Set == "Arctic") # Arctic An_Cont An_Pen

# Loop iterates over Cercozoa, diatoms, and green algae data
for(dataset in c("Cercozoa", "Diatoms", "GreenAlgae")){

  
  # ❄️ Load count table, metadata and taxonomy table
  counts <- read.table(paste0("CountData_", dataset, ".csv"), header = TRUE, sep = ";", dec = ".", row.names = 1)
  counts <- counts %>% select(ends_with("Sv")) # To select the Arctic
  # counts <- counts %>% select(starts_with("BIO")) # To select Continental An
  # counts <- counts %>% select(starts_with("C"), starts_with("K"), starts_with("Me")) # To select An Peninsula

  taxonomy <- read.table(paste0("Taxonomy_", dataset, ".csv"), header = TRUE, sep = ";", dec = ".", row.names = 1)
  
```  
## 2. Sparcity filtering
To reduce noise in the data, I addressed sparsity by filtering out OTUs that were present in too few samples. I grouped these rare OTUs into a single "pseudo taxon" to manage the sparsity issue, ensuring that only meaningful data were retained for further analysis. This step was crucial in refining the dataset and making it more robust for subsequent analyses.
``` r  
  print(dataset)
  # ❄️ Define a threshold for the prevalence filter
  threshold <- 5 # 10%: Arctic = 5; An_Pen = 5, An_Cont = 2, used 3 to filter out taxa without taxonomy 
  
  # Find all taxa that shall be grouped into one pseudeo-taxon
  to_filter <- which(rowSums(counts != 0) <= threshold)
  print(paste("OTUs present in <=", threshold, "of", ncol(counts),
              "samples are summed into one pseudo taxon."))
  
  # Filter counts data
  print(paste(length(to_filter)-nrow(counts[rowSums(counts) == 0,]), "out of",
              nrow(counts[rowSums(counts) != 0,]), "OTUs are filtered out."))
  print(paste0(sum(counts[to_filter,]), " out of ", sum(counts[-to_filter,])," (", 
               round(sum(counts[to_filter,])/sum(counts)*100,2),
               "%) reads are binned into the pseudo taxon."))
 
  counts_filtered <- counts[-to_filter,]
  counts <- rbind(counts_filtered, t(data.frame(filtered_taxa=colSums(counts[to_filter,]))))

  # Filter taxonomy data, save them in list
  taxonomy <- taxonomy[rownames(taxonomy) %in% rownames(counts),]
  TaxonomyData[[dataset]] <- taxonomy
```  

## 3. Compositionality data transformation 
I transformed the count data for compositional analysis by performing a centered log-ratio (clr) transformation. Before applying the transformation, I adjusted the data to replace zeros with a small constant to avoid computational issues. This transformation was essential for normalizing the data, allowing for accurate comparison across samples.  
``` r
# Transpose counts to have OTUs as columns
counts_transposed <- data.frame(t(counts))

# Replace 0 by 1/10 of 1
counts_transposed[counts_transposed == 0] <- 0.1

# clr transformation
counts_transposed_clr <- data.frame(clr(counts_transposed))

# Save counts in list
CountData[[dataset]] <- counts_transposed_clr

# Append transformed data to a list
Transposed_CountData[[dataset]] <- counts_transposed
}

 
 # Merge count data and taxonomy data
counts <- do.call(cbind, CountData)
taxonomy <- do.call(rbind, TaxonomyData)

Transposed_Counts_An_Cont <- do.call(cbind, Transposed_CountData)
write.csv(Transposed_Counts_An_Cont, file = "Per_Set_Thr10per100/CountsTransposed_ForNetwork_An_Cont.csv")

metadata$ID <- row.names(metadata)
metadata <-  dcast(metadata, ID ~ Set, length)
metadata <- metadata[,c("Arctic")] # ❄️ "Arctic", "An_Cont", "An_Pen", 

 

 # ❄️ Transform and export the count data 
  write.table(counts,
              file = "Per_Set_Thr10per100/FlashWeaveNetwork_Poles_CountData_Arctic.csv", sep=",", # Arctic An_Cont An_Pen
              quote=FALSE, row.names=F)
 
 # Export the metadata except Locations and Sampling_date 
 write.table(metadata,
             file = "FlashWeaveNetwork_Poles_Metadata.csv", sep=",", 
             quote=FALSE, row.names= F)  
 
 # ❄️ Export the taxonomy table
  write.csv(taxonomy, file = "Per_Set_Thr10per100/FlashWeaveNetwork_Poles_Taxonomy_Arctic.csv",  # Arctic An_Cont An_Pen
            quote=FALSE, row.names=T) 
``` r  
 
## 4. Run FlashWeave
 
FlashWeave is implemented in Julia, a high-level, high-performance dynamic language for technical computing. 

To run the script: 
1. Open a terminal or command prompt.
2. Start the Julia interpreter by running the following command:
 > /../julia-1.7.2/bin/julia
 
3. ❄️ Run the script: 
 > include("/../FlashWeave_NetworkInference_Arctic.jl")
 
The script is different for each region. The input data includes the count data filtered and transformed during the first three steps of this script. The following lines are written in a text file but saved as .jl in the used directory. Check the example available in [this directory](3_Bash)  

  
      ># Load FlashWeave
      > using FlashWeave
     
      ># Define input path for the count data
      > count_data_path = string("/home/alle/WORKING_DIR/Cristina/Networks/FlashWeaveNetwork_Poles_CountData.csv")
     
      ># Define input path for the metadata
      > meta_data_path = string("/home/alle/WORKING_DIR/Cristina/Networks/FlashWeaveNetwork_Poles_Metadata.csv")
     
      ># Calculate network
      > network = learn_network(count_data_path, meta_data_path, sensitive=true, heterogeneous=false)
     
      ># Save network
      > save_network(string("/home/alle/WORKING_DIR/Cristina/Networks/NetworkFlashWeave.edgelist"), network)
     
     
 
## 5. Transform FlashWeave network for Cytoscape            
  
### 5.1 Transform network data
``` r
# ❄️ Load data, change column names
 edges <- read.table('Per_Set_Thr10per100/FlashWeaveNetwork_Arctic.edgelist', sep="\t")  # Arctic.edgelist  An_Cont.edgelist  An_Pen.edgelist 
 colnames(edges) <- c("Node1", "Node2", "weight")
 
# Delete Cercozoa, Diatoms, or GreenAlgae from names
 edges$Node1 <- gsub("Cercozoa.|Diatoms.|GreenAlgae.", "", edges$Node1)
 edges$Node2 <- gsub("Cercozoa.|Diatoms.|GreenAlgae.", "", edges$Node2)
 
# Remove binned taxa from table (from prevalence filter)
edges <- edges[which(edges$Node1 != 'filtered_taxa' & edges$Node2 != 'filtered_taxa'),]
 
# Remove self loops
 edges <- edges[edges$Node1 != edges$Node2,]
 
 # Create new column specifying if correlation is positive or negative
 edges$Association <- "negative"
 edges$Association[edges$weight > 0] <- "positive"
 
 edges$Node1 <- sub("^X", "", edges$Node1)
 edges$Node2 <- sub("^X", "", edges$Node2)
 
 # Get number of negative and positive edges
 print(paste0("Number of positive edges: ", sum(edges$weight > 0)))
 print(paste0("Number of negative edges: ", sum(edges$weight < 0)))
 print(paste0("Total number of edges: ", nrow(edges)))
 
 # ❄️ Export data as .txt file
 write.table(edges, file = 'Per_Set_Thr10per100/Results/FlashWeaveNetwork_Arctic.txt', # Sp
             sep="\t", quote=FALSE, row.names=FALSE)
 ``` r
 

### 6. Create taxonomy table                    
``` r 
 # Create lists for taxonomy data
 TaxonomyData <- list()
 
 # Loop iterates over Cercozoa, diatoms, and green algae data
 for(dataset in c("Cercozoa", "Diatoms", "GreenAlgae")){
   
   # Load count table and taxonomy table
   counts <- read.table(paste0("CountData_", dataset, ".csv"), header = TRUE, sep = ";", dec = ".", row.names = 1)
   taxonomy <- read.table(paste0("Taxonomy_", dataset, ".csv"), header = TRUE, sep = ";", dec = ".", row.names = 1)
   
   # Merge number of reads per taxon with taxonomy table
   taxonomy <- merge(taxonomy, data.frame(NrReads=rowSums(counts)), by=0, all = T)
   
   # Rename rownames column to key
   colnames(taxonomy)[colnames(taxonomy) == "Row.names"] <- "Key"
   
   # Save counts and taxonomy in list
   TaxonomyData[[dataset]] <- taxonomy
 }
 
 # Merge taxonomy data
 taxonomy <- do.call(rbind, TaxonomyData)
 
 taxonomy$Key <- gsub("-", ".", taxonomy$Key)
 
 # Replace NAs
 taxonomy[is.na(taxonomy)] <- ""
 
 # Filter taxonomy, keep only taxa present in network
 taxonomy <- taxonomy[taxonomy$Key %in% c(edges$Node1, edges$Node2),]
 
 unique_nodes <- unique(c(edges$Node1, edges$Node2)) # 135 unique taxa, but I have only 133 in the taxonomy... why? 
 
 # ❄️ Export taxonomy data as txt file
 write.table(taxonomy, file = "Per_Set_Thr10per100/Results/FlashWeaveNetwork_Arctic_Taxonomy.txt",  # Sp
             sep="\t", quote=FALSE, row.names=FALSE)
``` 
 

### 7. Summarize networks
With this step I aggregated the networks for ease of visualization. The two final files are loded into Cytoscape for visualization. 
``` r
 
 # Install packages
 # if (!require("docstring")) install.packages("docstring")
 # if (!require("dplyr")) install.packages("dplyr")
 
 # Load packages
 library(docstring)
 library(dplyr)
 
 # Load association.frequencies function
 source("Functions/association.frequencies.R")
 
 # To see the documentation of this function we can use the package docstring 
 # docstring(association.frequencies)
 
 # Load edges and taxonomy table 
 edges <- read.table("Per_Set_Thr10per100/Results/FlashWeaveNetwork_Arctic.txt",
                     header = T, sep = "\t", dec = ".")
 taxonomy <- read.table("Per_Set_Thr10per100/Results/FlashWeaveNetwork_Arctic_Taxonomy.txt",
                        header = T, sep = "\t")
 
 entries_not_present <- setdiff(unique_nodes, taxonomy$Key)
 
 if (length(entries_not_present) == 0) {
   print("All entries in unique_nodes are present in taxonomy$Key.")
 } else {
   print("Entries not present in taxonomy$Key:")
   print(entries_not_present)
 }
 
 
 # Reduce the edges table to only the two node columns a  nd the edge weight
 # (this is required by the function 'association.frequencies')
 edges <- edges[,c("Node1", "Node2", "weight")]
 
 # Calculate frequencies per Genus
 summary <- association.frequencies(edge.table=edges, taxonomy=taxonomy,
                                    old.key="Key", new.key="Order")
 
 # ❄️ Export network and taxonomy data as txt file
 write.table(summary$Edges, file = 'Per_Set_Thr10per100/Results/FlashWeaveNetwork_Arctic_Order.txt',
             sep="\t", quote=FALSE, row.names=FALSE)
 write.table(summary$Taxonomy, file = "Per_Set_Thr10per100/Results/FlashWeaveNetwork_Arctic_Order_Taxonomy.txt",  
             sep="\t", quote=FALSE, row.names=FALSE)
```  