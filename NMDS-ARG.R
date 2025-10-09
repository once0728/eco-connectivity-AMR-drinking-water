# ---NMDS Visualisation ---
# Author: once0728
# Description: This script analyses differences in resistance gene composition between distinct sample groups.


# Set working directory
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# Load required packages
library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# Read ARG subtype abundance data
type <- read.table("ARG_subtype_abundance.txt", header = TRUE, row.names = 1, 
                   check.names = FALSE, sep = "\t", fileEncoding = 'GBK')

# Data preprocessing: remove NA values
type <- type %>%
  drop_na() %>%
  select_if(~ !all(is.na(.)))

# Transpose data matrix (samples as rows, genes as columns)
type <- t(type)

# Calculate Bray-Curtis dissimilarity matrix
type.distance <- vegdist(type, method = 'bray')

# Perform NMDS ordination analysis
set.seed(123) # Set random seed for reproducibility
df_nmds <- metaMDS(type.distance, k = 2, trymax = 100)
df_nmds_stress <- round(df_nmds$stress, 3)

# Extract NMDS coordinates
df_points <- as.data.frame(df_nmds$points)
df_points$samples <- rownames(df_points)
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')

# Read sample group information
group <- read.table("group.txt", sep = '\t', header = TRUE)
colnames(group) <- c("samples", "group")

# Merge NMDS coordinates with group information
df <- merge(df_points, group, by = "samples")

# Set group factor levels
df$group <- factor(df$group, levels = c("INF", "AS", "EFF", "RU", "RS", "RD", "SW", "TW"))

# Perform PERMANOVA to test significance of group differences
set.seed(123)
permanova <- adonis2(type.distance ~ group, data = df, permutations = 999)
permanova_r2 <- round(permanova$R2[1], 3)
permanova_p <- round(permanova$`Pr(>F)`[1], 3)

# Create color 
color <- c("#aa3a4b","#3b729b","#508870","#e29e4b","#cb7fab","#c0c05a","#60C1C6","#B69F89")

# Create visualization
p <- ggplot(data = df, aes(x = NMDS1, y = NMDS2, fill = group)) +
  theme_bw(base_size = 12) +
  geom_point(shape = 21, color = 'black', size = 3.5, alpha = 0.8) +
  
  # Add confidence ellipses (70% confidence interval)
  stat_ellipse(geom = "polygon", level = 0.70, 
               linetype = 2, size = 0.6, alpha = 0.2, aes(fill = group)) +
  
  # Add reference lines
  geom_vline(xintercept = 0, lty = "dashed", linewidth = 0.5, color = 'grey60') +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.5, color = 'grey60') +
  
  # Color settings
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  
  # Theme adjustments
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.8),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 11)
  ) +
  
  # Add statistical annotations
  labs(
    title = "Non-metric Multidimensional Scaling of ARG Profiles",
    subtitle = paste0("Stress = ", df_nmds_stress, 
                      "; PERMANOVA R² = ", permanova_r2, 
                      ", p = ", ifelse(permanova_p < 0.001, "< 0.001", permanova_p)),
    x = "NMDS1",
    y = "NMDS2"
  )

# Display plot
print(p)

# Save high-resolution figures (suitable for journal submission)
ggsave("NMDS_ARG.pdf", p, 
       width = 8, height = 6, dpi = 600)
ggsave("NMDS_ARG.png", p, 
       width = 8, height = 6, dpi = 600)