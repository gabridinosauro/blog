### Code for blog post Unkwnown


#BACTERIA/ARCHAEA ALPHA DIVERSITY
ASVtab = readRDS( "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Sky_islands_project/Tables/ASVtab_clean.RDS") 
TAXtab = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Sky_islands_project/Tables//TAXtab_clean.RDS")
metadata_bac = readRDS("//Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Sky_islands_project/Tables/metadata_bac.RDS")

####alpha diversity calculations
sort(colSums(ASVtab))
rar.ASVtab = rrarefy(t(ASVtab), 37000) #rarefaction 
sort(colSums(t(rar.ASVtab))) ##rarefaction worked
source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")
TAXtab
taxo = generate.tax.summary.modified(t(ASVtab), TAXtab)
Phylum_tab  = taxo$tax3
Class_tab = taxo$tax4
Order_tab  = taxo$tax5
Family_tab  = taxo$tax6
genus_tab  = taxo$tax7
genus_tab_rel = decostand(Family_tab, 2,method = "total")
colSums(genus_tab_rel) # ok
genus_tab_t = data.frame(t(genus_tab_rel))
rownames(metadata_bac) == rownames(genus_tab_t) # ok
metadata_bac$unknown = genus_tab_t$unknown
plot(metadata_bac$unknown~metadata_bac$Elevation)
cor.test(metadata_bac$unknown,metadata_bac$Elevation, method = "spearman")
richness.lmer<-lmerTest::lmer(unknown~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)






### Count how many singleton ASVs are present in each point
asv_counts <- ifelse(t_data > 0,1,0)
asv_counts = asv_counts[rowSums(asv_counts) ==1,]
metadata_bac$singletons = colSums(asv_counts)
plot(metadata_bac$singletons~metadata_bac$Elevation)
cor.test(metadata_bac$singletons,metadata_bac$Elevation, method = "spearman")
richness.lmer<-lmerTest::lmer(singletons~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)
r2(richness.lmer)
