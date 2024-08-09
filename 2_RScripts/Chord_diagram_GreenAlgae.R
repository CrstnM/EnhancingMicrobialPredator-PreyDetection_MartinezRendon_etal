# Chord diagram

# Cristina Martínez Rendón
# 14-06-2023


### Chord diagrams, as representations of presence/absence of genera in three geographical locations. Metabarcoding data.
### In the first 209 lines I process the data. You can download the "GreenA_chord.csv" file and jump directly to making the chord diagram and see how the chordDiagram() function works. Some steps might be redundant in those first lines.


library(tidyverse)
library(circlize)

rm(list = ls())

setwd("~/R_Projects/ArcticAntarctica/GreenAlgae") # Location of this R script. s

## 1. Data handling 
 
### Raw count table: table with OTUs as rows and individual samples as columns. The last column provides the taxonomic assignments for each OTU.
counttable <- read.table("305WP2PolarGA.unique.agc.txt", header=T,row.names = 1) 

    # Relativized count table. Non-rarefied data, filtered (K2|TOT|TCtr|C1 sites), deleted OTUs that col.summed=0, or those present in =< 3 samples. I need both tables, this one does not have the taxonomy column. This table I produced in the NMDS script.
    relcounttable <- read.table("Species_mat_relativ.txt", header=T,row.names = 1) 

      
### Filter out things from the raw table using base R. 
counttable <- counttable[,-which(colnames(counttable) %in% c("SUM","repSeqName","repSeq","mock_community"))]
counttable <- as.data.frame(t(counttable))
row.names(counttable) <- gsub("^X","",row.names(counttable))
counttable <- filter(counttable, !grepl( "K2|TOT|TCtr|unused|C1", row.names(counttable))) 

counttable <- as_tibble(counttable, rownames = "names") %>%
  mutate(names_site = str_replace(names, "(\\d_\\d)(?!.*_An)", "\\1_An")) %>% 
  select(names_site, starts_with("Otu"))
counttable <- as.data.frame(counttable)
row.names(counttable) <- counttable[,1]
counttable <- counttable[, -1]
counttablet <- t(counttable)
counttablet <- as.data.frame(counttablet)

        
    ### Add the taxonomy column to the relativized count table
    relcounttablet <- t(relcounttable)    # rarecounttablet <- t(rarecounttable)
    relcounttablet <- as.data.frame(relcounttablet)
    
    common_rows <- data.frame(rownames(counttablet) %in% rownames(relcounttablet))
    subset <- counttablet[unlist(common_rows), ]
    relcounttablet$OTUConTaxonomy <- subset$OTUConTaxonomy
    
    
# Split the strings in OTUConTaxonomy, separated by ";"
taxonomy_split <- as_tibble(relcounttablet, rownames="OTU") %>% 
  select(OTU, OTUConTaxonomy) %>% 
  separate(OTUConTaxonomy, sep=";", into = c("Kingdom", "Supergroup", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  mutate_all(~str_replace_all(., "\\(100\\)", ""))
taxonomy_split <- as.data.frame(taxonomy_split)


# I wanted to present the results in phylogenetic order, so I sorted the data:
taxonomy_split$Order <- factor(taxonomy_split$Order, levels = c("Monomastigales", "Chlorellales", "Prasiolales", "Microthamniales",
                                                                  "Trebouxiales", "Watanabeales", "Trebouxiophyceae", "Trentepohliales", "Ulotrichales",
                                                                  "Ulvales", "Scotinosphaerales", "Chaetophorales", "Sphaeropleales", "Chlamydomonadales"))
taxonomy_split$Class <- factor(taxonomy_split$Class, levels = c("Mamiellophyceae", "Trebouxiophyceae", "Ulvophyceae", "Chlorophyceae"))
taxonomy_split <- arrange(taxonomy_split, Class, Order, Genus, Species)
taxonomy_split <- taxonomy_split %>% 
  mutate_all(~str_replace_all(., "_X{1,2}", "") %>% 
               str_replace_all("-", "_"))

### Here I downloaded the taxonomy_split and kept on editing in excel. I entered exact family names, many were falsely assigned! 
  # This can often happen with the reference databases. Edited spp. names when needed (Unk_x)
# write.table(taxonomy_split, "taxonomy_split", sep = "\t", quote = FALSE, row.names = TRUE)

# Reload the taxonomy_split object. Reorganize it by OTU number here, so later on the wrong taxonomy won't be assigned to a OTU.
taxonomy_split <- read.delim("taxonomy_split2.txt", header=T,row.names = 1)
          
taxonomy_split <- arrange(taxonomy_split, ON) %>% 
  select(-ON)

# 2. Generate a sorted table with the unique taxa
taxonomy_unique <- taxonomy_split %>% 
  select(-OTU, -Species) %>% 
  distinct()

taxonomy_unique$Class <- factor(taxonomy_unique$Class, levels = c("Mamiellophyceae", "Trebouxiophyceae", "Chlorophyceae", "Ulvophyceae"))

taxonomy_unique$Order <- factor(taxonomy_unique$Order, levels = c("Monomastigales", "Chlorellales", "Prasiolales", "Microthamniales",
                                                                  "Trebouxiales", "Watanabeales", "Trebouxiophyceae", "Trentepohliales", "Ulotrichales",
                                                                  "Ulvales", "Scotinosphaerales", "Chaetophorales", "Sphaeropleales", "Chlamydomonadales"))

taxonomy_unique <- arrange(taxonomy_unique, Class, Order)


## 👌😉  ##
  

# 3. Split the three sets.

# Integrate columns from the metadata info and subset the samples
info <- read.table("SampleMetadata.txt", header = TRUE, row.names = 1)
relcounttablent <- relcounttable[row.names(relcounttable) != "OTUConTaxonomy", ]

    common_rows <- data.frame(rownames(info) %in% rownames(relcounttablent))
    subset <- info[unlist(common_rows), ]
    relcounttablent$set <- subset$set
    
relcounttablent <- relcounttablent[, c("set", names(relcounttablent)[-ncol(relcounttablent)])]

# Subset for the three regions while deleting those OTU lines that are not present in each set. 
Arctic_sub <- relcounttablent %>% 
  filter(set == "Arctic") %>% 
  select(-set) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(as.numeric) %>% 
  mutate(Totals = rowSums(.)) %>% 
  mutate(Arctic = Totals > 0) %>% 
  select(Arctic)
Arctic_sub$Class <- as.factor(taxonomy_split$Class)
Arctic_sub$Order <- as.factor(taxonomy_split$Order)
Arctic_sub$Family <- as.factor(taxonomy_split$Family)
Arctic_sub$Genus <- taxonomy_split$Genus
Arctic_sub <- Arctic_sub %>% 
  filter(Arctic == TRUE)
  
An_Pen_sub <- relcounttablent %>% 
  filter(set == "An_Pen") %>% 
  select(-set) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(as.numeric) %>% 
  mutate(Totals = rowSums(.)) %>% 
  mutate(An_Pen = Totals > 0) %>% 
  select(An_Pen)
An_Pen_sub$Class <- as.factor(taxonomy_split$Class)
An_Pen_sub$Order <- as.factor(taxonomy_split$Order)
An_Pen_sub$Family <- as.factor(taxonomy_split$Family)
An_Pen_sub$Genus <- taxonomy_split$Genus
An_Pen_sub <- An_Pen_sub %>% 
  filter(An_Pen == TRUE)

An_Cont_sub <- relcounttablent %>% 
  filter(set == "An_Cont") %>% 
  select(-set) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(as.numeric) %>% 
  mutate(Totals = rowSums(.)) %>% 
  mutate(An_Cont = Totals > 0) %>% 
  select(An_Cont)
An_Cont_sub$Class <- as.factor(taxonomy_split$Class)
An_Cont_sub$Order <- as.factor(taxonomy_split$Order)
An_Cont_sub$Family <- as.factor(taxonomy_split$Family)
An_Cont_sub$Genus <- taxonomy_split$Genus
An_Cont_sub <- An_Cont_sub %>% 
  filter(An_Cont == TRUE)

  
# 4. Fill the taxonomy_unique table with the occurrence of each unique OTU in each set
# Create a new column in taxonomy_unique
taxonomy_unique$Arctic <- FALSE
taxonomy_unique$An_Pen <- FALSE
taxonomy_unique$An_Cont <- FALSE


# Check if the four columns match with the values in Arctic_sub
for (i in 1:nrow(taxonomy_unique)) {
  for (j in 1:nrow(Arctic_sub)) {
    if (all(taxonomy_unique[i, c("Order", "Family", "Genus")] == Arctic_sub[j, c("Order", "Family", "Genus")])) {
      taxonomy_unique[i, "Arctic"] <- TRUE
      break  # exit inner loop once a match is found
    }
  }
}

for (i in 1:nrow(taxonomy_unique)) {
  for (j in 1:nrow(An_Pen_sub)) {
    if (all(taxonomy_unique[i, c("Order", "Family", "Genus")] == An_Pen_sub[j, c("Order", "Family", "Genus")])) {
      taxonomy_unique[i, "An_Pen"] <- TRUE
      break  # exit inner loop once a match is found
    }
  }
}

for (i in 1:nrow(taxonomy_unique)) {
  for (j in 1:nrow(An_Cont_sub)) {
    if (all(taxonomy_unique[i, c("Order", "Family", "Genus")] == An_Cont_sub[j, c("Order", "Family", "Genus")])) {
      taxonomy_unique[i, "An_Cont"] <- TRUE
      break  # exit inner loop once a match is found
    }
  }
}


# Add the colors each Class will have in the chord!
taxonomy_us_edited <- taxonomy_unique %>%
  mutate(Colors = case_when(
    Class == "Mamiellophyceae" ~ "#CCEBC5",
    Class == "Trebouxiophyceae" ~ "#78c679",
    Class == "Chlorophyceae" ~ "#35978F",
    Class == "Ulvophyceae" ~ "#016c59",
    TRUE ~ NA_character_
  ))

GreenA_chord <- taxonomy_us_edited %>% 
  dplyr::select(Genus, Arctic, An_Pen, An_Cont) %>% 
  mutate_all(~str_replace_all(., "TRUE", "1")) %>% 
  mutate_all(~str_replace_all(., "FALSE", "0")) %>% 
  mutate_at(vars(c("Arctic", "An_Cont", "An_Pen")), as.numeric)
  
row.names(GreenA_chord) <- GreenA_chord[,1]
GreenA_chord <- GreenA_chord[, -1]
GreenA_chord <- as.matrix(GreenA_chord) # Circlize uses a matrix as input.

      # I made this table available, if you wanted to jump directly to this line and avoid the processing. 
      # write.table(GreenA_chord, "GreenA_chord.csv", sep = "\t", quote = FALSE, row.names = TRUE)
      GreenA_chord <- as.matrix(read.table("GreenA_chord.csv", header=T,row.names = 1))

      
## 2. PLOT IT!
### Define colors
polar_colors <- c("Arctic" = "#FFE4B5", "An_Cont" = "#4A80BD", "An_Pen" = "#87B8D9")
taxa_colors <- setNames(taxonomy_us_edited$Color, taxonomy_us_edited$Genus)
combined_colors <- c(taxa_colors, polar_colors)


### Create the chord diagram without labels or axes
chordDiagram(GreenA_chord, annotationTrack = "grid", preAllocateTracks = 1, grid.col = combined_colors, 
             directional = 1, transparency = 0.2)

### Add the labels and not the axis ticks and numbers
circos.trackPlotRegion(track.index = 2, panel.fun = function(x,y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
# Print sector labels
  circos.text(mean(xlim), ylim[1] + 2.5, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex=0.6)
}, bg.border = NA)

### If I want to print the axis ticks and numbers... but I don't
# circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, 
#             sector.index = sector.name, track.index = 2)
# }, bg.border = NA)


### Save 
dev.copy(pdf, "Chord_withCirclize_GreenA5.pdf", width = 10, height = 9) # png, units="in", res=500)
dev.off()


circos.clear()




