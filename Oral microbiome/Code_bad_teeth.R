### Blog post analysis my oral microbiome
library(ggplot2)
library(ggpubr)

setwd("/Users/gabri/Library/CloudStorage/Box-Box/github_code/blog/Oral microbiome")
metaphlan_results = read.delim("LS149_S16_profiled_metagenome.txt", header = FALSE,
                               col.names = c("Taxonomy", "NCBI_tax_id", 
                                             "relative_abundance","additional_species"))
metaphlan_results = metaphlan_results[-c(1:5),]
metaphlan_results$Subject = "Gabri"
metaphlan_results$relative_abundance = as.numeric(metaphlan_results$relative_abundance)

mypalette2  <-c("#40004b","#ffffbf","#762a83","#de77ae","#f46d43","#f7f7f7","#d9f0d3","#a6dba0","#a6cee3","#1f78b4","#b2df8a","#33a02c","#4d9221",
                "#2166ac","#5aae61","#1b7837","#d73027","#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#b2182b","#de77ae","#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30","#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99")

### Select Phylum
metaphlan_results_phylum = metaphlan_results[grep("p__", metaphlan_results$Taxonomy),]
metaphlan_results_phylum = metaphlan_results_phylum[-grep("c__", metaphlan_results_phylum$Taxonomy),]
metaphlan_results_phylum$Taxonomy = sub(".*\\|p__(.*)", "\\1", metaphlan_results_phylum$Taxonomy)
a = ggplot(metaphlan_results_phylum, aes(fill=Taxonomy, y=relative_abundance, x= Subject)) + 
  geom_bar( stat="identity") +
  theme_bw() + scale_fill_manual(values = mypalette2) + labs(fill='Phylum') + 
  ylab("Relative Abundance") +xlab("Sample") 
a


### Select Species 
metaphlan_results_species = metaphlan_results[grep("s__", metaphlan_results$Taxonomy),]
metaphlan_results_species$Taxonomy = sub(".*\\|s__(.*)", "\\1", metaphlan_results_species$Taxonomy)
metaphlan_results_species = metaphlan_results_species[-grep("t__", metaphlan_results_species$Taxonomy),]
metaphlan_results_species50 = metaphlan_results_species[1:50,]
metaphlan_results_species50$Taxonomy = factor(metaphlan_results_species50$Taxonomy, levels = metaphlan_results_species50$Taxonomy)

b = ggplot(metaphlan_results_species50, aes(y=relative_abundance, x= Taxonomy)) + 
  geom_bar( stat="identity") +
  theme_bw() + scale_fill_manual(values = mypalette2) + labs(fill='Phylum') + 
  ylab("Relative Abundance") +xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

b



### Select Candidatus_Absconditabacteria
Candidatus_Absconditabacteria = metaphlan_results[grep("Candidatus_Absconditabacteria", metaphlan_results$Taxonomy),]
metaphlan_results_species = Candidatus_Absconditabacteria[grep("s__", Candidatus_Absconditabacteria$Taxonomy),]
metaphlan_results_species$Taxonomy = sub(".*\\|s__(.*)", "\\1", metaphlan_results_species$Taxonomy)
metaphlan_results_species = metaphlan_results_species[-grep("t__", metaphlan_results_species$Taxonomy),]
metaphlan_results_species50 = metaphlan_results_species#[1:50,]
metaphlan_results_species50$Taxonomy = factor(metaphlan_results_species50$Taxonomy, levels = metaphlan_results_species50$Taxonomy)

b = ggplot(metaphlan_results_species50, aes(y=relative_abundance, x= Taxonomy)) + 
  geom_bar( stat="identity") +
  theme_bw() + scale_fill_manual(values = mypalette2) + labs(fill='Phylum') + 
  ylab("Relative Abundance (%)") +xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Candidatus_Absconditabacteria")

b
Saccharibacteria
### Select Saccharibacteria
Saccharibacteria = metaphlan_results[grep("Saccharibacteria", metaphlan_results$Taxonomy),]
metaphlan_results_species = Saccharibacteria[grep("s__", Saccharibacteria$Taxonomy),]
metaphlan_results_species$Taxonomy = sub(".*\\|s__(.*)", "\\1", metaphlan_results_species$Taxonomy)
metaphlan_results_species = metaphlan_results_species[-grep("t__", metaphlan_results_species$Taxonomy),]
metaphlan_results_species50 = metaphlan_results_species#[1:50,]
metaphlan_results_species50$Taxonomy = factor(metaphlan_results_species50$Taxonomy, levels = metaphlan_results_species50$Taxonomy)

b = ggplot(metaphlan_results_species50, aes(y=relative_abundance, x= Taxonomy)) + 
  geom_bar( stat="identity") +
  theme_bw() + scale_fill_manual(values = mypalette2) + labs(fill='Phylum') + 
  ylab("Relative Abundance (%)") +xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Saccharibacteria")

b
