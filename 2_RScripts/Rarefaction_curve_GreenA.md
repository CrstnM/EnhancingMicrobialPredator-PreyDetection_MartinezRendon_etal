# Rarefaction curves for each sample


Cristina Martínez Rendón
17.10.2023

rm(list=ls())

**R version:** 4.3.0 (2023-04-21)

**Packages**

``` r
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggthemes")) install.packages("ggthemes")
if (!require("vegan")) install.packages("vegan")

library(tidyverse)
library(ggthemes) 
library(vegan)
``` r

setwd("~/R_Projects/ArcticAntarctica/GreenAlgae")

counts <- read.delim("305WP2PolarGA.unique.agc.txt", header=T,row.names = 1)
counts <- counts[,-which(colnames(counts) %in% c("SUM","repSeqName","repSeq","OTUConTaxonomy", "mock_community"))]
counts <- as.data.frame(t(counts))
counts <- filter(counts, !grepl( "K2|TOT|TCtr|unused", row.names(counts)))
row.names(counts) <- gsub("^X","",row.names(counts))

#Continue with a tibble
counts <- as_tibble(counts, rownames = "names") %>%
  mutate(names_site = str_replace(names, "(\\d_\\d)(?!.*_An)", "\\1_An")) %>% 
  select(names_site, everything()) %>% 
  select(-names)
counts <- as.data.frame(counts)
rownames(counts) <- counts$names_site
counts <- counts[,-1]


# Calculate rarefaction curves
# rarecurves <- rarecurve(counts, step = 50)
rarecurves_Cerco <- rarecurve(counts, step = 50, xlab = "Sample Size", ylab = "Species", label = FALSE) # Base R rarecurve.

# This line already creates the data frame usable for ggplot.
rare_counts <- rarecurve(counts, step = 50, tidy=TRUE)

# Load metadata
metadata <- read.delim("SampleMetadata.txt", header = TRUE, row.names = 1)

# Combine rarefaction data and metadata to one data frame to color-code by WWTP compartment 
# For by.y choose the column with your sample IDs, if by.y=0 rownames are used
# rarefaction_data <- merge(rarefaction_data, metadata, by.x="SampleID", by.y = 0, all = T)
rare_counts <- merge(rare_counts, metadata, by.x="Site", by.y = 0, all = T)



color <- c("#FFE4B5", "#225895", "#89aec2")

# Plot
  ggplot(rare_counts, aes(x = Sample, y = Species, group = Site, color = set)) +
  geom_line(size = 0.8) +
  # scale_color_manual(values = c(polar_colors)) +
  scale_color_manual(name = "Sampling sites",
                       values = c("#FFE4B5", "#225895", "#89aec2"),
                       breaks = c("Arctic", "An_Cont", "An_Pen")) +
  # Set plot title, x-axis label & y-axis label 
  labs(title = "Rarefaction curves per samplig sites, Green Algae", x = "Number of reads", y = "Number of OTUs") +
  # Change background of plot 
  # geom_vline(xintercept = 4800, color="gray40", size=1) +
  theme_minimal() +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12, face = "bold"), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.position = "right",
        legend.text.align = 0,
        legend.direction = "vertical", 
        strip.text.y = element_text(size=12, face = "bold"), 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size=12,face="bold")) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

ggsave(file = "Plots/Rarefaction_GreenAlgae.png", dpi=300, width = 8, height = 5)