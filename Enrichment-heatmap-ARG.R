# --- Heatmap Visualization of Antibiotic Resistance Gene Distribution and Enrichment ---
# Author: once0728
# Description: Generates a heatmap showing ARG distribution and enrichment across water types
# (wastewater, river, drinking water). Highlights enriched ARG types compared to baseline samples.

# =========================================================
# 1. Load Required Packages
# =========================================================
packages <- c("tidyverse", "RColorBrewer", "rstudioapi")
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

# Read and preprocess data
df <- read.csv("Total_ARG.csv", header = TRUE, stringsAsFactors = FALSE, 
               check.names = FALSE, fileEncoding = 'GBK')

# Transform to long format for ggplot compatibility
df_long <- pivot_longer(df,
                        cols = -c("type"),
                        names_to = "sample",
                        values_to = "value")

# Add sample grouping information
df_long$group1 <- case_when(
  grepl("INF|AS|EFF", df_long$sample) ~ "Wastewater",
  grepl("RU|RS|RD", df_long$sample) ~ "River water",
  grepl("SW|TW", df_long$sample) ~ "Drinking water"
)
df_long$group1 <- factor(df_long$group1, levels = c("Wastewater", "River water", "Drinking water"))

# Extract sample type prefix
df_long$sample_type <- substr(df_long$sample, 1, 2)

# Calculate baseline values for enrichment analysis
inf_base <- df_long %>% 
  filter(sample_type == "INF") %>%
  group_by(type) %>%
  summarise(INF_base = mean(value))

ru_base <- df_long %>% 
  filter(sample_type == "RU") %>%
  group_by(type) %>%
  summarise(RU_base = mean(value))

sw_base <- df_long %>% 
  filter(sample_type == "SW") %>%
  group_by(type) %>%
  summarise(SW_base = mean(value))

# Merge baseline values and calculate enrichment factors
df_long1 <- df_long %>%
  left_join(inf_base, by = "type") %>%
  left_join(ru_base, by = "type") %>%
  left_join(sw_base, by = "type") %>%
  mutate(
    enrichment = case_when(
      sample_type == "EFF" ~ (value - INF_base)/INF_base,
      sample_type == "RD" ~ (value - RU_base)/RU_base,
      sample_type == "TW" ~ (value - SW_base)/SW_base,
      TRUE ~ NA_real_
    ),
    symbol = case_when(
      enrichment > 1 & enrichment <= 10 ~ "+",
      enrichment > 10 ~ "++",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-INF_base, -RU_base, -SW_base)

# Set factor levels for proper ordering
df_long1$sample <- factor(df_long1$sample, levels = c("EFF", "AS", "INF", "RD", "RS", "RU", "TW", "SW"))
df_long1$type <- factor(df_long1$type, levels = c("Beta-lactam", "Bacitracin", "MLS", "Aminoglycoside", "Sulfonamide", 
                                                  "Tetracycline", "Polymyxin", "Chloramphenicol", "Trimethoprim", "Quinolone", 
                                                  "Rifamycin", "Florfenicol", "Other peptide antibiotics", "Fosfomycin", 
                                                  "Mupirocin", "Bicyclomycin", "Novobiocin", "Streptothricin", "Bleomycin", 
                                                  "Vancomycin", "Pleuromutilin tiamulin", "Antibacterial fatty acid", "Puromycin", 
                                                  "Defensin", "Edeine"))

# Replace symbol codes with appropriate characters
df_long1[which(df_long1$symbol == "E"), "symbol"] <- '+'   
df_long1[which(df_long1$symbol == "EE"), "symbol"] <- '++'

# Define color palette (diverging blue-white-red)
heatmap_colors <- c("#5493c4", "white", "#cd5f5e")

# Create heatmap visualization
p <- ggplot(df_long1, aes(type, sample)) +
  geom_tile(aes(fill = log10(value + 0.000001)), color = "#E8E8E8", size = 0.5) +
  scale_fill_gradientn(
    colours = heatmap_colors,
    name = expression(paste("ARG abundance (log"[10], "(copy/cell))")),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 12,
      barheight = 0.8,
      direction = "horizontal",
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  facet_grid(group1 ~ ., scales = "free", space = 'free', switch = "y") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label = symbol), 
            vjust = 0.8,
            size = 4,
            fontface = "bold") +
  labs(
    x = " ",
    y = " ")+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 15)),
    plot.caption = element_text(size = 9, hjust = 0, face = "italic"),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# Display plot
print(p)

# Save high-resolution figures
ggsave("Figure_Heatmap_ARG_Distribution.pdf", p, width = 10, height = 7)
ggsave("Figure_Heatmap_ARG_Distribution.png", p, width = 10, height = 7)

