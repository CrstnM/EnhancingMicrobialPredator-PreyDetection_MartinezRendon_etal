# Alpha diversity: Diversity indexes using vegan 

# Cristina Martínez Rendón
# 13-06-2023


library(tidyverse) check
library(vegan)     check
library(cowplot)
library(iNEXT)
library(agricolae)
library(FSA)
library(multcompView)
library(car)

rm(list = ls())

setwd("~/R_Projects/ArcticAntarctica/GreenAlgae") # Location of this R script. 

set.seed(850511)

# 1. Data handling

# Non-rarefied data, filtered (K2|TOT|TCtr|C1 sites), deleted OTUs that col.summed=0, or those with =< 3 
counttable <- read.table("305WP2PolarGA.unique.agc.txt", header=T,row.names = 1)

    # Relativized count table
      relcounttable <- read.table("Species_mat_relativ.txt", header=T,row.names = 1) 

    counttable <- relcounttable
          

####### - HILL NUMBERS - #######
 # 1. OTU Richness: This corresponds to 𝑞= 0, where no weighting by relative abundance occurs, so every OTU is counted equally, representing the total species (or OTU) richness.
    
 # 2. Exponential Shannon: This corresponds to 𝑞= 1, which is the exponential of the Shannon entropy. This metric takes into account both the richness and evenness of species in the community, giving more weight to species that are more evenly distributed.
    
 # 3. Inverse Simpson: This corresponds to 𝑞= 2, which is the inverse of the Simpson index. This metric gives more weight to the more abundant species, reflecting the probability that two randomly chosen individuals from the dataset belong to different species.
    
    
# We calculate Hill numbers with alpha = 0, 1 & 2 from count data for each sample using the function renyi (package vegan).
# If needed, we transpose the count data to have sample IDs as the row names and taxa IDs as the column names.

# Calculate hill numbers for alpha = 0, 1, 2
renyi <- renyi(counttable, hill = T, scales = c(0, 1, 2)) # I ran this with the relativized count table. I tried it with the count table with
# non-rarefied abundances, with almost the exact same result as to the created numbers. 
renyi$Set <- factor(metadata_subset$set, levels = c("Arctic", "An_Pen", "An_Cont"))
renyi$uniqsample <- as.factor(metadata_subset$uniqsample)

# Is each level normally distributed?
hill0 <- shapiro.test(renyi$"0")
hill1 <- shapiro.test(renyi$"1")
hill2 <- shapiro.test(renyi$"2")

# Print the test results for each index column
print(hill0)     # p-value = 0.0008843 ... Kruskal-Wallis, then Wilcoxon's and Dunn's tests
print(hill1)     # p-value = 0.09894 ... Meets normality! ANOVA? Heterogeneous variances, then Wilcoxon's and Dunn's tests
print(hill2)     # p-value = 0.02279 ... Kruskal-Wallis, then Wilcoxon's and Dunn's tests

# Change the column names in the renyi object
colnames(renyi) <- c("hill0", "hill1", "hill2", "Set", "uniqsample")

          ### First test the homogeneity of variances
          leveneTest(hill1 ~ Set, data = renyi)     # Pr 0.001233 **
          # The variances show no homogeneity. Use alternative methods that do not assume equal variances, such as non-parametric tests.


##### NON PARAMETRIC TESTS #####
          
# Kruskal-Wallis
kruskal.test(hill0  ~ Set, data = renyi)
# data:  otu_richness by Set
# Kruskal-Wallis chi-squared = 18.097, df = 2, p-value = 0.2539
            hill0_groups <- data.frame(
            groups = c("a", "a", "a"),
            Set = factor(c("An_Cont", "Arctic", "An_Pen")))
kruskal.test(hill1 ~ Set, data = renyi)
# data:  shannon by Set
# Kruskal-Wallis chi-squared = 15.051, df = 2, p-value = 0.001179 * 
kruskal.test(hill2 ~ Set, data = renyi)
# data:  simpson by Set
# Kruskal-Wallis chi-squared = 8.9269, df = 2, p-value = 0.01976 *


# Post-hoc pairwise comparison:

      # Hill 1

dunn_hill1 <- dunnTest(hill1 ~ Set, data = renyi, method = "bonferroni")
        #       Comparison         Z      P.unadj        P.adj
        # 1 An_Cont - An_Pen -3.6145964 0.0003008159 0.0009024476   * different
        # 2 An_Cont - Arctic -3.0089396 0.0026216126 0.0078648378   * different
        # 3  An_Pen - Arctic  0.7090686 0.4782819215 1.0000000000   * similar

wilcox_hill1 <- pairwise.wilcox.test(renyi$hill1, renyi$Set, p.adjust.method = "bonferroni")
      # Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
      #         Arctic An_Pen
      # An_Pen  0.8028 -     
      # An_Cont 0.0014 0.0013
      # P value adjustment method: bonferroni 

# To express the p-values as letters:

#### For dunn_hill1

# Extract the adjusted p-values
p_adj_hill1 <- dunn_hill1$res$P.adj
names(p_adj_hill1) <- dunn_hill1$res$Comparison
      # > p_adj_hill1
      # An_Cont - An_Pen An_Cont - Arctic  An_Pen - Arctic 
      # 0.0009024476     0.0078648378     1.0000000000 

# Apply multcompLetters() to get letter groups
letter_groups_hill1 <- multcompLetters(p_adj_hill1)
      # letter_groups_hill1                Not too bad, but now I have two letters for An_Pen... (again)
      # An_Cont   An_Pen    An_Pen   Arctic 
      # "a"     "ab"      "b"      "b" 


#### For wilcox_hill1

p_adj_hill1 <- wilcox_hill1$p.value # Extract the adjusted p-values
      # > p_adj_hill1              
      # Arctic      An_Pen
      # An_Pen  1.000000000           NA
      # An_Cont 0.007262721 0.0009146086

# Convert lower triangular matrix to long format
melted_whill1 <- p_adj_hill1 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Comparison") %>%
  pivot_longer(cols = -Comparison, names_to = "Comparison2", values_to = "P_Value") %>%
  filter(!is.na(P_Value)) %>% 
  unite(Comparison, Comparison, Comparison2, sep = " - ") # This merges the first two columns

# Convert melted to a named vector
named_vector_whill1 <- deframe(melted_whill1)

# Apply multcompLetters() to get letter groups
letter_groups_whill1 <- multcompLetters(named_vector_whill1)
      # > letter_groups_whill1 ................. same result as with dunn_shannon
      # An_Pen  An_Cont    Arctic   An_Pen 
      # "ab"      "a"      "b"      "b" 

# Create hill1_groups to plot 
hill1_groups <- data.frame(groups = letter_groups_whill1$Letters) # Extract the letters from the list object
hill1_groups$Set <- rownames(hill1_groups)   # Create a new column with the Sets, otherwise the object is degraded to a vector in the next step
hill1_groups <- hill1_groups[-1, ]           # Remove duplicate entry for "An_Pen"
hill1_groups$Set <- trimws(hill1_groups$Set) # Remove leading/trailing white space in names
hill1_groups$Set <- factor(hill1_groups$Set, levels = c("Arctic", "An_Pen", "An_Cont"))
hill1_groups$rgroups <- c("b", "a", "a")



# Hill 2
dunn_hill2 <- dunnTest(hill2 ~ Set, data = renyi, method = "bonferroni")
      #       Comparison         Z      P.unadj        P.adj
      # 1 An_Cont - An_Pen -2.7694263 0.00561551 0.01684653   * different
      # 2 An_Cont - Arctic -2.2522481 0.02430660 0.07291979   * similar
      # 3  An_Pen - Arctic  0.6130944 0.53981392 1.00000000   * similar

wilcox_hill2 <- pairwise.wilcox.test(renyi$hill2, renyi$Set, p.adjust.method = "bonferroni")
      # Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
      #         Arctic An_Pen
      # An_Pen  1.000  -     
      # An_Cont 0.070  0.018 
      # P value adjustment method: bonferroni 

# To express the p-values as letters:

#### For dunn_hill2

# Extract the adjusted p-values
p_adj_hill2 <- dunn_hill2$res$P.adj
names(p_adj_hill2) <- dunn_hill2$res$Comparison
      # > p_adj_hill2
      # An_Cont - An_Pen An_Cont - Arctic  An_Pen - Arctic 
      # 0.01684653       0.07291979       1.00000000 

# Apply multcompLetters() to get letter groups
letter_groups_hill2 <- multcompLetters(p_adj_hill2)
      # letter_groups_hill2               
      # An_Cont   An_Pen    An_Pen   Arctic 
      # "a"     "ab"      "b"      "ab" 


#### For wilcox_hill2

p_adj_hill2 <- wilcox_hill2$p.value # Extract the adjusted p-values
      # > p_adj_hill2              
      # Arctic      An_Pen
      # An_Pen  1.00000000         NA
      # An_Cont 0.07017737 0.01832549

# Convert lower triangular matrix to long format
melted_whill2 <- p_adj_hill2 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Comparison") %>%
  pivot_longer(cols = -Comparison, names_to = "Comparison2", values_to = "P_Value") %>%
  filter(!is.na(P_Value)) %>% 
  unite(Comparison, Comparison, Comparison2, sep = " - ") # This merges the first two columns

# Convert melted to a named vector
named_vector_whill2 <- deframe(melted_whill2)

# Apply multcompLetters() to get letter groups
letter_groups_whill2 <- multcompLetters(named_vector_whill2)
      # > letter_groups_whill2 ................. same result as with dunn_shannon
      # An_Pen  An_Cont    Arctic   An_Pen 
      # "ab"      "a"      "ab"      "b" 

# Create hill2_groups to plot 
hill2_groups <- data.frame(groups = letter_groups_whill2$Letters) # Extract the letters from the list object
hill2_groups$Set <- rownames(hill2_groups)   # Create a new column with the Sets, otherwise the object is degraded to a vector in the next step
hill2_groups <- hill2_groups[-1, ]           # Remove duplicate entry for "An_Pen"
hill2_groups$Set <- trimws(hill2_groups$Set) # Remove leading/trailing white space in names
hill2_groups$Set <- factor(hill2_groups$Set, levels = c("Arctic", "An_Pen", "An_Cont"))
# hill2_groups$rgroups <- c("ab", "b", "a")



color <- c("#E1C697", "#5F93B5", "#225895") 
color_light <- c("#FFE4B5", "#4A80BD", "#87B8D9")


hill0_plot <- ggplot(renyi, aes(x = Set, y = hill0, fill = Set)) +
  geom_boxplot(width = 0.8, color = color, size = 0.8) +
  scale_fill_manual(values = color_light,
                    limits = c("Arctic", "An_Cont", "An_Pen")) +
  geom_jitter(aes(colour = Set), width = 0.05, height = 0.02, size = 2) +
  scale_color_manual(values = color) +  # Apply the color values to the jitter layer
  labs(y = "Hill 0") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12, face = "bold"), 
        axis.title.x=element_blank(),
        legend.position = "none", # originally, "right". Here I deleted the legend to plot the three alpha indexes together.
        strip.text.y = element_text(size=12, face = "bold")) +
  theme(strip.text.x = element_text(size=12,face="bold")) +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_text(data = hill0_groups,
            aes(x = Set, y = Inf, label = groups), vjust = 1, size = 5); # Add group labels 
hill0_plot


hill1_plot <- ggplot(renyi, aes(x = Set, y = hill1, fill = Set)) +
  geom_boxplot(width = 0.8, color = color, size = 0.8) +
  scale_fill_manual(values = color_light,
                    limits = c("Arctic", "An_Cont", "An_Pen")) +
  geom_jitter(aes(colour = Set), width = 0.05, height = 0.02, size = 2) +
  scale_color_manual(values = color) +  # Apply the color values to the jitter layer
  labs(y = "Hill 1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12, face = "bold"), 
        axis.title.x=element_blank(),
        legend.position = "none", # originally, "right". Here I deleted the legend to plot the three alpha indexes together.
        strip.text.y = element_text(size=12, face = "bold")) +
  theme(strip.text.x = element_text(size=12,face="bold")) +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_text(data = shannon_groups,
            aes(x = Set, y = Inf, label = rgroups), vjust = 1, size = 5); # Add group labels 
hill1_plot


hill2_plot <- ggplot(renyi, aes(x = Set, y = hill2, fill = Set)) +
  geom_boxplot(width = 0.8, color = color, size = 0.8) +
  scale_fill_manual(values = color_light,
                    limits = c("Arctic", "An_Cont", "An_Pen")) +
  geom_jitter(aes(colour = Set), width = 0.05, height = 0.02, size = 2) +
  scale_color_manual(values = color) +  # Apply the color values to the jitter layer
  labs(y = "Hill 2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12, face = "bold"), 
        axis.title.x=element_blank(),
        legend.position = "none", # originally, "right". Here I deleted the legend to plot the three alpha indexes together.
        strip.text.y = element_text(size=12, face = "bold")) +
  theme(strip.text.x = element_text(size=12,face="bold")) +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_text(data = shannon_groups,
            aes(x = Set, y = Inf, label = rgroups), vjust = 1, size = 5); # Add group labels 
hill2_plot

plot_grid(hill0_plot, hill1_plot, hill2_plot, ncol=3)

ggsave(file = "Plots/Alpha_diversity_Diatoms_HillNumbers_116samples_relativ.png", dpi=300, width = 15, height = 5)


