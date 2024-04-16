#BACTERIA/ARCHAEA ALPHA DIVERSITY
library(ggplot2)
library(phyloseq)

ASVtab = readRDS( "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Sky_islands_project/Tables/ASVtab_clean.RDS") 
TAXtab = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Sky_islands_project/Tables//TAXtab_clean.RDS")
metadata_bac = readRDS("//Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Sky_islands_project/Tables/metadata_bac.RDS")
# you can find the same files in the github repository



####alpha diversity calculations
sort(colSums(ASVtab))
rar.ASVtab = rrarefy(t(ASVtab), 37000) #rarefaction 
sort(colSums(t(rar.ASVtab))) ##rarefaction worked
source("https://raw.githubusercontent.com/gabridinosauro/Useful_functions/main/helpful_r_functions.R")


taxo = generate.tax.summary.modified(t(ASVtab), TAXtab)
Phylum_tab  = taxo$tax3
tab_rel = decostand(Phylum_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(metadata_bac) == rownames(tab_rel_t) # ok
metadata_bac$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)
phyl = ggplot(metadata_bac, aes(x = Elevation, y = log(unknown), color = Mount )) +
  geom_point() +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("Phylum - unknown")

Class_tab  = taxo$tax4
tab_rel = decostand(Class_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(metadata_bac) == rownames(tab_rel_t) # ok
metadata_bac$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)
class = ggplot(metadata_bac, aes(x = Elevation, y = log(unknown), color = Mount )) +
  geom_point() +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("Class - unknown")
class

order_tab  = taxo$tax5
tab_rel = decostand(order_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(metadata_bac) == rownames(tab_rel_t) # ok
metadata_bac$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)
order = ggplot(metadata_bac, aes(x = Elevation, y = log(unknown), color = Mount )) +
  geom_point() +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("Order - unknown")
order

family_tab  = taxo$tax6
tab_rel = decostand(family_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(metadata_bac) == rownames(tab_rel_t) # ok
metadata_bac$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)
fam = ggplot(metadata_bac, aes(x = Elevation, y = log(unknown), color = Mount )) +
  geom_point() +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("Family - unknown")
fam


genus_tab  = taxo$tax7
tab_rel = decostand(genus_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(metadata_bac) == rownames(tab_rel_t) # ok
metadata_bac$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Elevation + (1|Mount), data = metadata_bac)
summary(richness.lmer)
anova(richness.lmer)
gen = ggplot(metadata_bac, aes(x = Elevation, y = log(unknown), color = Mount )) +
  geom_point() +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("Genus - unknown")
gen



ggarrange(phyl, class, order, fam, gen, common.legend = TRUE)






### Dust project ------
`%nin%` = negate(`%in%`)
ASVtab_dust<-readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/phyloseq_dust16sASVnames.RDS")
ASVtab_dust=prune_samples(sample_data(ASVtab_dust)$point!=26, ASVtab_dust)
#ASVtab_dust <- prune_taxa(taxa_sums(ASVtab_dust) > 0, ASVtab_dust)
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine SS"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))
asvtab_dust = data.frame(ASVtab_dust@otu_table)
asvtab_dust = data.frame(t(asvtab_dust))
taxa = data.frame(ASVtab_dust@tax_table)
sam_data_dust = sam_data_dust[which(sam_data_dust$Dust_Type %nin% c("Coarse SS", "Medium SS", "Fine SS")),]
rownames(sam_data_dust) = gsub("-",".",rownames(sam_data_dust))
asvtab_dust = asvtab_dust[,colnames(asvtab_dust) %in% rownames(sam_data_dust)]
asvtab_dust= asvtab_dust[,match(rownames(sam_data_dust),colnames(asvtab_dust))]
rownames(sam_data_dust)== colnames(asvtab_dust) # ok
taxo = generate.tax.summary.modified(t(asvtab_dust), taxa)

taxo = generate.tax.summary.modified(t(ASVtab), TAXtab)
Phylum_tab  = taxo$tax3
tab_rel = decostand(Phylum_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(sam_data_dust) == rownames(tab_rel_t) # ok
sam_data_dust$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Dust_Type + (1|Point_type), data = sam_data_dust)
summary(richness.lmer)
anova(richness.lmer)
Phylum = ggplot(sam_data_dust, aes(x = Dust_Type, y = unknown)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Point_type )) +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
   theme_bw() + ggtitle("Phylum - unknown")
Phylum


Class_tab  = taxo$tax4
tab_rel = decostand(Class_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(sam_data_dust) == rownames(tab_rel_t) # ok
sam_data_dust$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Dust_Type + (1|Point_type), data = sam_data_dust)
summary(richness.lmer)
anova(richness.lmer)
Class = ggplot(sam_data_dust, aes(x = Dust_Type, y = unknown)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Point_type )) +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  theme_bw() + ggtitle("Class - unknown")
Class



Order_tab  = taxo$tax5
tab_rel = decostand(Order_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(sam_data_dust) == rownames(tab_rel_t) # ok
sam_data_dust$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Dust_Type + (1|Point_type), data = sam_data_dust)
summary(richness.lmer)
anova(richness.lmer)
Order = ggplot(sam_data_dust, aes(x = Dust_Type, y = unknown)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Point_type )) +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  theme_bw() + ggtitle("Order - unknown")
Order


Family_tab  = taxo$tax6
tab_rel = decostand(Family_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(sam_data_dust) == rownames(tab_rel_t) # ok
sam_data_dust$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Dust_Type + (1|Point_type), data = sam_data_dust)
summary(richness.lmer)
anova(richness.lmer)
Family = ggplot(sam_data_dust, aes(x = Dust_Type, y = unknown)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Point_type )) +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  theme_bw() + ggtitle("Family - unknown")
Family



Genus_tab  = taxo$tax7
tab_rel = decostand(Genus_tab, 2,method = "total")
colSums(tab_rel) # ok
tab_rel_t = data.frame(t(tab_rel))
rownames(sam_data_dust) == rownames(tab_rel_t) # ok
sam_data_dust$unknown = tab_rel_t$unknown
richness.lmer<-lmerTest::lmer(log(unknown)~Dust_Type + (1|Point_type), data = sam_data_dust)
summary(richness.lmer)
anova(richness.lmer)
Genus = ggplot(sam_data_dust, aes(x = Dust_Type, y = unknown)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Point_type )) +
  labs(x = "Elevation", y = "log(Relative abundance)") + 
  theme_bw() + ggtitle("Genus - unknown")
Genus



ggarrange(Phylum, Class, Order, Family, Genus, common.legend = TRUE)
