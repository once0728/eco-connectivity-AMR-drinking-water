# --- Antibiotic Resistome Risk Score Analysis ---
# Author: once0728
# Description: This script analyzes the relationship between antibiotic resistance gene (ARG) 
# abundance and resistome risk scores.

# Load packages
packages <- c("ggplot2", "cowplot", "tidyverse", "ggsci", "ggpmisc", "readxl", "ggpubr")
for(p in packages){
  if(!require(p, character.only = TRUE)){
    install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
}

# Set working directory
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# Read and prepare data
df <- read_excel("RiskScore-ARG_MGE.xlsx")
head(df)

# Filter for ARG data only
df1 <- df %>% filter(type1 == "ARG")


# Define color palette
col <- c("#aa3a4b", "#3b729b")

# Create scatter plot with regression lines
p <- ggplot(df1, aes(x = x, y = y, fill = type2)) +
  geom_point(
    shape = 21,
    size = 4,
    alpha = 0.7,
    color = "white",
    stroke = 0.5
  ) +
  geom_smooth(
    method = "lm",
    aes(color = type2),
    se = TRUE,
    formula = y ~ x,
    linetype = 1,
    alpha = 0.2,
    size = 1
  ) +
  stat_poly_eq(
    formula = y ~ x,
    size = 4,
    aes(
      color = type2,
      label = paste(..rr.label.., ..p.value.label.., sep = "~~~")
    ),
    parse = TRUE,
    label.x = "right",
    label.y = "top"
  ) +
  facet_wrap(
    . ~ type2,
    nrow = 2,
    ncol = 2,
    scales = "free"
  ) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.title.x = element_text(colour = "black", size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(colour = "black", size = 13, face = "bold", margin = margin(r = 10)),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
    plot.caption = element_text(size = 9, hjust = 0, face = "italic")
  ) 


print(p)


# Save high-resolution figures
ggsave("Figure_Resistome_Risk_Scores_ARG.pdf", p, width = 10, height = 8, dpi = 600)
ggsave("Figure_Resistome_Risk_Scores_ARG.png", p, width = 10, height = 8, dpi = 600)