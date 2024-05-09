### PERMANOVA results simulation
library(ggplot2)
library(tidyr)
library(ecole)
### Create dummy distance matrix
# Set the seed for reproducibility
set.seed(123)
### create a dummy mapping file
mapping_file = data.frame(sample = c(1:120), group = c(rep("Before",40), rep("Treatment_A",40), rep("Treatment_B",40)), pair = c(1:40, 1:40, 1:40))


#### Create a distance matrix
distance_matrix <- matrix(0, nrow = 120, ncol = 120)
# Assign distances between each of the first 50 and the latter 50
for (i in 1:120) {
  for (j in 1:120) {
    distance_matrix[i, j] <- runif(1, 0, 1)
    # Random number between 0 and 1
  }
}

### reduce distances within groups (make them cluster)
distance_matrix[1:40, 1:40] = distance_matrix[1:40, 1:40] - 0.5*distance_matrix[1:40, 1:40]
distance_matrix[41:80, 41:80] = distance_matrix[41:80, 41:80] - 0.8*distance_matrix[41:80, 41:80]
distance_matrix[81:120, 81:120] = distance_matrix[81:120, 81:120] - 0.5*distance_matrix[81:120, 81:120]

### Increase distance between before and B
for(i in 1:40) 
{
  distance_matrix[80+i, i] = distance_matrix[80+i, i] +  2*distance_matrix[80+i, i]
}
distance_matrix = ifelse(distance_matrix > 1, 1, distance_matrix)



bac.nmds<-metaMDS(distance_matrix, k=2, try=100)
bac.nmds$stress #0.06
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

mapping_file$Axis01<-bac.nmds$points[,1]
mapping_file$Axis02<-bac.nmds$points[,2]

### mean distance between pairs
bac_nmds1<-ggplot(mapping_file, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=group), size = 3)+
  theme_classic() +
  geom_line(aes(group = pair), alpha = 0.1)  
bac_nmds1

#mean distance between pairs
mapping_file_dist = mapping_file[1:40,]
mapping_file_dist$Scenario_1 = NA
for(i in 1:40) 
{
  value = distance_matrix[40+i, i]
  mapping_file_dist$Treatment_A[i] = value
  value = distance_matrix[80+i, i]
  mapping_file_dist$Treatment_B[i] = value
}


mapping_file_long <- pivot_longer(mapping_file_dist, 
                                  cols = c(Treatment_A, Treatment_B), 
                                  names_to = "distribution", 
                                  values_to = "value")



d = ggplot(mapping_file_long, aes(x = distribution, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distances between groups", 
       x = "Distribution", 
       y = "BC - dissimilarities before vs treatment")

d

disper = betadisper(as.dist(distance_matrix),mapping_file$group)
mapping_file$distances = disper$distances
e = ggplot(mapping_file, aes(x = group, y = distances)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distances within groups", 
       x = "group", 
       y = "BC - distances from centroids")

e



ggarrange(bac_nmds1, d, e, widths = c(2,1,1), ncol = 3, common.legend = TRUE)
permanova.pairwise = permanova_pairwise(distance_matrix, c(mapping_file$group),  permutations = 999,distance = "bray", padj = "fdr")
permanova.pairwise


#mean distance between pairs
mapping_file_dist = mapping_file[1:40,]
mapping_file_dist$Scenario_1 = NA
for(i in 1:50) 
{
  value = distance_matrix[50+i, i]
  mapping_file_dist$Scenario_1[i] = value
}



perm = adonis(distance_matrix ~group, data =mapping_file , strata = mapping_file$pair )
perm$aov.tab


bac.nmds<-metaMDS(distance_matrix, k=2, try=100)
bac.nmds$stress #0.06
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

mapping_file$Axis01<-bac.nmds$points[,1]
mapping_file$Axis02<-bac.nmds$points[,2]

### mean distance between pairs
bac_nmds1<-ggplot(mapping_file, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=group), size = 3)+
  theme_classic()+
  annotate("text",x=min(mapping_file$Axis01)+ 0.2,y=max(mapping_file$Axis02),hjust=1,label= "PERMANOVA r2 = 0.01")+
  geom_line(aes(group = pair), alpha = 0.1)  + ggtitle("Scenario 1")
bac_nmds1

