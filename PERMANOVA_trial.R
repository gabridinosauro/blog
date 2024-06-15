### PERMANOVA results simulation
library(ggplot2)
library(tidyr)
### Create dummy distance matrix
# Set the seed for reproducibility
set.seed(123)
mapping_file = data.frame(sample = c(1:100), group = c(rep("A",50), rep("B",50)), pair = c(1:50, 1:50))

# Number of samples
n <- 100

### try with paired samples
distance_matrix <- matrix(0, nrow = n, ncol = n)
# Assign distances between each of the first 50 and the latter 50
for (i in 1:100) {
  for (j in 1:100) {
    distance_matrix[i, j] <- runif(1, 0, 1)
    # Random number between 0 and 1
  }
}


#mean distance between pairs
mapping_file_dist = mapping_file[1:50,]
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




### Scenario2 - increase distances between paired points
distance_matrix1 = distance_matrix
distance_matrix2 = as.dist(distance_matrix1)
distance_matrix1 = as.matrix(distance_matrix2)

for(i in 1:50) 
{
  distance_matrix1[50+i, i] = distance_matrix1[50+i, i] +  2*distance_matrix1[50+i, i]
}
distance_matrix1 = ifelse(distance_matrix1 > 1, 1, distance_matrix1)

perm = adonis(distance_matrix1 ~group, data =mapping_file, strata = mapping_file$pair )
perm$aov.tab


bac.nmds<-metaMDS(distance_matrix1, k=2, try=100)
bac.nmds$stress #0.06
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

mapping_file$Axis01<-bac.nmds$points[,1]
mapping_file$Axis02<-bac.nmds$points[,2]



bac_nmds2<-ggplot(mapping_file, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=group), size = 3)+
  theme_classic()+
  geom_line(aes(group = pair), alpha = 0.1)+
  annotate("text",x=min(mapping_file$Axis01)+ 0.2,y=max(mapping_file$Axis02),hjust=1,label= "PERMANOVA r2 = 0.009")+
  ggtitle("Scenario 2")
bac_nmds2




mapping_file_dist$Scenario_2 = NA
for(i in 1:50) 
{
  value = distance_matrix1[50+i, i]
  mapping_file_dist$Scenario_2[i] = value
}

mapping_file_long <- pivot_longer(mapping_file_dist, 
                                  cols = c(Scenario_1, Scenario_2), 
                                  names_to = "distribution", 
                                  values_to = "value")


ggplot(mapping_file_long, aes(x = distribution, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Distances", 
       x = "Distribution", 
       y = "Value")












###
distance_matrix2 = distance_matrix
distance_matrix2[1:50, 1:50] = distance_matrix2[1:50, 1:50] - 0.5*distance_matrix2[1:50, 1:50]
distance_matrix2[51:100, 51:100] = distance_matrix2[51:100, 51:100] - 0.5*distance_matrix2[51:100, 51:100]
perm = adonis(distance_matrix2 ~group, data =mapping_file, strata = mapping_file$pair )
perm$aov.tab


bac.nmds<-metaMDS(distance_matrix2, k=2, try=100)
bac.nmds$stress #0.06
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

mapping_file$Axis01<-bac.nmds$points[,1]
mapping_file$Axis02<-bac.nmds$points[,2]



bac_nmds3<-ggplot(mapping_file, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=group), size = 3)+
  theme_classic()+
  geom_line(aes(group = pair), alpha = 0.1)+
  annotate("text",x=min(mapping_file$Axis01)+ 0.2,y=max(mapping_file$Axis02),hjust=1,label= "PERMANOVA r2 = 0.27***") +
  ggtitle("Scenario 3")
bac_nmds3




mapping_file_dist$Scenario_3 = NA
for(i in 1:50) 
{
  value = distance_matrix2[50+i, i]
  mapping_file_dist$Scenario_3[i] = value
}

mapping_file_long <- pivot_longer(mapping_file_dist, 
                                  cols = c(Scenario_1, Scenario_2, Scenario_3), 
                                  names_to = "distribution", 
                                  values_to = "value")
d = ggplot(mapping_file_long, aes(x = distribution, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Distances", 
       x = "Distribution", 
       y = "BC - dissimilarities")


ggarrange(bac_nmds1,bac_nmds2,bac_nmds3,d, common.legend = TRUE)






