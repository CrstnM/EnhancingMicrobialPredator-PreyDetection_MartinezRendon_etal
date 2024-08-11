Enhancing Microbial Predator-Prey Detection with Network and Trait-Based Analyses
================

This is a collection of the scripts detailing the methods of metabarcoding data analyses used for the paper in process authored by Cristina Martínez Rendón, Christina Braun, Maria Kappelsberger, Jens Boy, Angélica Casanova-Katny, Karin Glaser, and Kenneth Dumack.

After the release of the publication, the raw sequencing data will be available to download in [GenBank](https://www.ncbi.nlm.nih.gov/sra/PRJNA1144814), under the accession numbers SAMN43044305, SAMN43044306, and SAMN43044307 for Cercozoa, green algae, and ochrophytes, respectively. The three original OTU tables with taxonomic annotation are available in the folder [1_Data](1_Data), where you can also access the metadata, including chemical edaphic parameters. Intermediate files and results can be obtained by running the scripts, while the figures produced are given at the end of each script. Here, we present the scripts applied for the green algal data set, which were applied identically to the other two datasets.

It should be noted that some of the figure parameters may vary somewhat after reanalysis owing to the random iterative nature of certain analyses (such as rarefaction, beta diversity, and PCoA), but the fundamental conclusions will largely stay the same.


Sequence processing
----------
We produced independent metabarcoding datasets of a predator group (Cercozoa) and their respective, putative prey in polar biocrusts (here green algae and ochrophytes). First, after generating the OTU tables, we [rarefied](2_RScripts/Rarefaction_curve_GreenAlgae.md) the data to check for saturation, which was achieved for the three datasets. Downstream analyses were then conducted with the non-rarefied datasets.

Diversity measures
----------
Community composition was visualized through [chord plots](2_RScripts/Chord_diagram_GreenAlgae.md) and [Venn diagrams](2_RScripts/Venn_diagram_GreenAlgae.md). [Alpha](2_RScripts/Alpha_GreenAlgae.md) and [beta diversity](2_RScripts/NMDS_GreenAlgae.md) metrics were as well visualized and statistically tested. 

Cross-kingdom co-occurrence network inference
----------
Biotic putative interactions were detected and visualized via the calculation of [cross-kingdom co-occurrence networks](2_RScripts/Cross-kingdom_networks) among the three investigated taxa. To reduce spurious putative interactions, and based on the NMDS analyses visualized from beta diversity dissimilarities, we conducted separate network analyses for each sampling region. Moreover, prevalence filters were applied to exclude rare taxa. Finally, significant putative interactions between taxa (nodes) were aggregated at the order level. The input data (raw data without taxonomy) and metadata for the network analyses are available [here](https://github.com/CrstnM/EnhancingMicrobialPredator-PreyDetection_MartinezRendon_etal/tree/main/1_Data/Network_Analyses]
