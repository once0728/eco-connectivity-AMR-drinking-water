# --- Bubble Plot Visualization of Antibiotic-Resistant and Multidrug-Resistant Bacterial Strains ---
# Author: once0728
# Description: Generates a bubble plot visualization showing the proportion of antibiotic-resistant 
# and multidrug-resistant bacterial strains across different sample types.
# Date: 2025-10-09

# --------------------- Environment Setup ---------------------
packages <- c("ggplot2", "tidyverse", "ggpubr", "RColorBrewer", "rstudioapi")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Set working directory automatically (if in RStudio)
if ("rstudioapi" %in% installed.packages()) {
  dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(dir)
}


# Set working directory
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# Read and preprocess data
data <- read.csv("DATA.csv", sep = ",", stringsAsFactors = FALSE,
                 check.names = FALSE, fileEncoding = 'GBK')


# Transform data to long format
df_long <- pivot_longer(data,
                        cols = -c("ID"),
                        names_to = "Sample",
                        values_to = "value")

# Define color palette
col <- c("#aa3a4b", "#3b729b", "#508870", "#e29e4b", "#cb7fab", "#c0c05a", "#60C1C6", "#B69F89")

# Set factor levels for proper ordering
df_long$Sample <- factor(df_long$Sample, levels = c("INF", "AS", "EFF", "RU", "RS", "RD", "SW", "TW"))
df_long$ID <- factor(df_long$ID, levels = rev(c("Amp", "Tet", "Kan", "Smz", "AMP-KAN", "AMP-TET", "AMP-SMZ", 
                                                "KAN-TET", "KAN-SMZ", "TET-SMZ", "AMP-KAN-TET", "AMP-KAN-SMZ", 
                                                "AMP-TET-SMZ", "KAN-TET-SMZ", "AMP-KAN-TET-SMZ")))

# Create bubble plot
p2 <- ggplot(df_long, aes(x = Sample, y = ID, size = value, color = Sample)) +
  # Plot non-zero values as solid circles
  geom_point(
    data = subset(df_long, value != 0),
    alpha = 0.85,
    shape = 16
  ) +
  # Plot zero values as hollow circles
  geom_point(
    data = subset(df_long, value == 0),
    alpha = 0.85,
    size = 3,
    shape = 21,
    fill = NA,
    color = "grey70"
  ) +
  # Scale settings
  scale_size(
    range = c(1, 10),
    name = "Proportion of\nResistant Strains",
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
    labels = c("0.1", "0.3", "0.5", "0.7", "0.9")
  ) +
  scale_color_manual(
    values = col,
    name = "Sample Type"
  ) +
  # Theme and styling
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linetype = "dashed", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      colour = "black",
      size = 10
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 10
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold",
      margin = margin(r = 10)
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      face = "italic"
    )
  ) +
  # Guides for legends
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    size = guide_legend()
  )

# Display plot
print(p2)

ggsave('Figure_Antibiotic_Resistance_Prevalence.pdf', p2, 
       width = 10, height = 7, dpi = 600)
ggsave('Figure_Antibiotic_Resistance_Prevalence.tiff', p2, 
       width = 10, height = 7, dpi = 600)

