library(ggplot2)
library(corrplot)


### code for post genomic ratios
MAGs_data = readRDS("MAGs_data.RDS")
KEGG_level1 = readRDS("counts_level1_MAGs.RDS")

# 1. Compute the correlation matrix
cor_matrix <- cor(KEGG_level1, method = "spearman", use = "pairwise.complete.obs")

# Create the heatmap
corrplot(cor_matrix, 
         method = "color", 
         type = "full",  # Show the full matrix
         tl.cex = 0.7,   # Reduce text size for labels
         tl.col = "black", 
         addCoef.col = NA,  # Remove correlation coefficients for simplicity
         cl.cex = 0.8)   # Adjust color legend text size


rownames(KEGG_level1) == MAGs_data$genome
KEGG_level1_rel = KEGG_level1/MAGs_data$totalORFs

cor_matrix <- cor(KEGG_level1_rel, method = "spearman", use = "pairwise.complete.obs")

# Create the heatmap
corrplot(cor_matrix, 
         method = "color", 
         type = "full",  # Show the full matrix
         tl.cex = 0.7,   # Reduce text size for labels
         tl.col = "black", 
         addCoef.col = NA,  # Remove correlation coefficients for simplicity
         cl.cex = 0.8)   # Adjust color legend text size


library(ggpubr)

# Create the plot with correlation line and p-value
a = ggplot(KEGG_level1_rel, aes(x = Genetic.Information.Processing, 
                            y = Environmental.Information.Processing)) +
  geom_point(color = "darkblue") +  # Add scatter points
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +  # Add Pearson regression line (no confidence interval shading)
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.15) +  # Add correlation coefficient and p-value
  theme_classic() +  # Clean theme
  labs(x = "fraction of ORFs annotated to Genetic Information Processing", 
       y = "fraction of ORFs annotated to Environmental Information Processing")

b = ggplot(KEGG_level1_rel, aes(x = Genetic.Information.Processing, 
                            y = Metabolism)) +
  geom_point(color = "darkblue") +  # Add scatter points
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +  # Add Pearson regression line (no confidence interval shading)
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.15) +  # Add correlation coefficient and p-value
  theme_classic() +  # Clean theme
  labs(x = "fraction of ORFs annotated to Genetic Information Processing", 
       y = "fraction of ORFs annotated to Metabolism")

ggarrange(a,b)



# Combine the data into a single data frame for ggplot
data_combined <- data.frame(
  sums = rep(MAGs_data$totalORFs, 3),  # Repeat sums for both datasets
  processing = c(KEGG_level1$Genetic.Information.Processing, 
                 KEGG_level1$Environmental.Information.Processing,
                 KEGG_level1$Metabolism),  # Combine processing data
  type = rep(c("Genetic Information Processing", "Environmental information processing", "Metabolism"), each = length(MAGs_data$totalORFs))  # Label each type
)

# Create the scatterplot
ggplot(data_combined, aes(x = sums, y = processing, color = type)) +
  geom_point(size = 1, alpha = 0.5) +  # Add points
  labs(
    x = "Total number of ORFs",
    y = "KO count",
    color = "Type"
  ) +
  theme_classic() +
  geom_smooth(method = "lm")# Use a clean theme

