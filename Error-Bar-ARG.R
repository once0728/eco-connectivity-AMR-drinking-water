# --- Bar Plot Visualization of Total Antibiotic Resistance Gene Abundance ---
# Author: once0728
# Revised: 2025-10-09
# Description: Generates a publication-quality bar plot showing total antibiotic resistance gene (ARG)
# abundance across sample types, with error bars and standardized color palette.

# =========================================================
# 1. Load Required Packages
# =========================================================
packages <- c("ggplot2", "dplyr", "scales")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}



# =========================================================
# 2. Working Directory Setup
# =========================================================
# Set working directory
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# Read and prepare data
df <- read.csv('Total_ARG.csv')
head(df)

# Set factor levels for proper ordering
df$ID <- factor(df$ID, levels = c("INF", "AS", "EFF", "RU", "RS", "RD", "SW", "TW"))

# Define color 
col <- c("#aa3a4b","#3b729b","#508870","#e29e4b","#cb7fab","#c0c05a","#60C1C6","#B69F89")

# Create bar plot with error bars
p <- ggplot(df, aes(x = ID)) +
  geom_bar(aes(y = Mean, fill = ID),  # Using fill instead of color for bars
           color = "black", size = 0.5, 
           stat = "identity", width = 0.7, alpha = 0.8) +  
  scale_fill_manual(values = col) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.2, size = 0.5, color = "black") +
  labs(y = "ARG Abundance", x = "Sample Groups",
       title = "Total Antibiotic Resistance Gene Abundance Across Sample Types") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.8),
    axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.y = element_text(size = 12, colour = "black", face = "bold"),
    axis.title.x = element_text(size = 12, colour = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "none"  
  )

# Display plot
print(p)

# Save high-resolution figures
ggsave('Total_ARG.pdf', p, width = 6, height = 5, dpi = 600)
ggsave('Total_ARG.png', p, width = 6, height = 5, dpi = 600)