#BACTERIA/ARCHAEA ALPHA DIVERSITY
library(ggplot2)


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
