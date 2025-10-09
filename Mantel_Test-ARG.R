# --- Mantel Test Analysis of ARGs and Environmental Factors ---
# Author: once0728
# Revised: 2025-10-09
# Description: Performs Mantel test between antibiotic resistance genes (ARGs)
# and environmental parameters, visualized as a network-enhanced correlation heatmap.

# =========================================================
# 1. Load Required Packages
# =========================================================
packages <- c("dplyr", "linkET", "ggplot2", "RColorBrewer", "cols4all")
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
df <- read.csv("arg_subtype_all-richness.csv", sep = ",", header = TRUE, 
               row.names = 1, check.names = FALSE)
env <- read.csv("ENV_all.csv", sep = ",", header = TRUE, 
                row.names = 1, check.names = FALSE)


# Perform Mantel test
mantel_result <- mantel_test(df, env,
                             spec_select = list(
                               "ARG" = 1:1309,
                               "MGE" = 1310:1351
                             ))

# Process Mantel test results for visualization
mantel_processed <- mantel_result %>% 
  mutate(
    r.sign = cut(r, breaks = c(-Inf, 0, Inf), 
                 labels = c("Negative", "Positive")),
    p.sign = cut(p, breaks = c(0, 0.05, Inf), 
                 labels = c("P < 0.05", "P ≥ 0.05"),
                 include.lowest = TRUE,
                 right = FALSE),
    r.abs = cut(abs(r), breaks = c(-Inf, 0.25, 0.5, Inf),
                labels = c("< 0.25", "0.25 - 0.5", "≥ 0.5"),
                include.lowest = TRUE,
                right = FALSE)
  )

# Display Mantel test results
cat("\nMantel test results:\n")
print(mantel_processed)

# Calculate environmental factor correlations
env_cor <- correlate(env, method = "spearman")

# Create base correlation heatmap
p_base <- qcorrplot(env_cor,
                    grid_col = "grey50",
                    grid_size = 0.2,
                    type = "upper",
                    diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    name = "Spearman's ρ"
  ) +
  theme(
    axis.text = element_text(color = "black", size = 10, face = "bold"),
    axis.text.x.top = element_text(
      color = "black", size = 10, angle = 45, 
      hjust = 0, vjust = 0, face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# Add significance markers
p_sig <- p_base +
  geom_mark(
    size = 3.5,
    only_mark = TRUE,
    sig_level = c(0.05, 0.01, 0.001),
    sig_thres = 0.05,
    mark = c("*", "**", "***"),
    colour = "black"
  )

# Add Mantel test connections
p_connected <- p_sig +
  geom_couple(
    data = mantel_processed,
    aes(
      colour = r.sign,
      size = r.abs,
      linetype = p.sign
    ),
    curvature = 0.2,
    nudge_x = 0.2,
    label.fontface = 2,
    label.family = "sans",
    label.size = 4
  )

# Final styling and legend customization
p_final <- p_connected +
  scale_size_manual(
    values = c("< 0.25" = 0.6, "0.25 - 0.5" = 1.0, "≥ 0.5" = 1.6),
    name = "Mantel |r|"
  ) +
  scale_colour_manual(
    values = c("Negative" = "#4a97c5", "Positive" = "#d0765c"),
    name = "Correlation\nDirection"
  ) +
  scale_linetype_manual(
    values = c("P < 0.05" = "solid", "P ≥ 0.05" = "dashed"),
    name = "Significance"
  ) +
  guides(
    fill = guide_colorbar(
      title = "Spearman's ρ",
      barwidth = 1,
      barheight = 8,
      order = 1
    ),
    linetype = guide_legend(
      title = "Mantel p-value",
      override.aes = list(size = 1),
      order = 3
    ),
    colour = guide_legend(
      title = "Correlation\nDirection",
      override.aes = list(size = 1),
      order = 4
    ),
    size = guide_legend(
      title = "Mantel |r|",
      override.aes = list(colour = "black"),
      order = 2
    )
  ) +
  theme(
    legend.key = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.text = element_text(color = "black", size = 9),
    legend.title = element_text(color = "black", size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 9, hjust = 0)
  ) 

# Display final plot
print(p_final)

# Save high-resolution figures
ggsave("Mantel_Test.pdf", p_final, width = 10, height = 8)
ggsave("Mantel_Test.tiff", p_final, width = 10, height = 8)

