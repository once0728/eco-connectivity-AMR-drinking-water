# --- Integrated Visualization of Antibiotic Resistance Genes in Metagenome-Assembled Genomes(MAGs) ---
# Author:  once0728
# Description: This script generates a comprehensive multi-panel visualization of antibiotic resistance genes (ARGs)
# and mobile genetic elements (MGEs) in MAGs. 


# Load required packages
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(patchwork)
library(RColorBrewer)

# Set working directory
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# Read and process phylogenetic tree
tree <- read.tree("bin.unrooted.tree")

# Read annotation data
anno <- readxl::read_xlsx("anno_arg.xlsx")

# Define color palette
mycol <- c("#c13736", "#2978b5", "#508870", "#e29e4b", "#cb7fab", "#c0c05a", "#60C1C6", "#B69F89")
mylevel <- c("Proteobacteria (23)", "Actinobacteriota (4)", "Bacteroidota (3)",
             "Firmicutes (2)", "Firmicutes_C (2)", "Chloroflexota (2)",
             "Firmicutes_A (1)", "Verrucomicrobiota (1)")

# Prepare phylum data for tree coloring
df_phylum <- subset(anno, select = c(ID, Phylum))
list_phylum <- split(df_phylum$ID, df_phylum$Phylum)
tree_phylum <- groupOTU(tree, list_phylum)

# Create phylogenetic tree visualization
p1 <- ggtree(tree_phylum, layout = "rectangular", alpha = 0.45, size = 0.8) +
  scale_color_manual(values = mycol, na.value = "#000000", guide = "none") +
  geom_tiplab(align = TRUE, size = 3.5) +
  theme(plot.margin = margin(5, 5, 5, 5))

# Add phylum annotation to the tree
p2 <- p1 + new_scale_fill() +
  geom_fruit(
    data = anno,
    geom = geom_tile,
    mapping = aes(y = ID, fill = factor(Phylum, levels = mylevel)),
    alpha = 0.7,
    pwidth = 0.06,
    offset = -1.04
  ) +
  scale_fill_manual(
    values = mycol,
    na.translate = FALSE,
    guide = guide_legend(
      keywidth = 0.8,
      keyheight = 0.8,
      ncol = 1,
      title = "Phylum",
      override.aes = list(color = "white", size = 4),
      order = 3
    )
  ) +
  theme(
    legend.position = c(0.25, 0.89),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# Read and process ARG abundance data for bubble plot
df <- read.csv("bin_arg.csv", header = TRUE, sep = ",")
df_long <- pivot_longer(
  df,
  cols = -c("ID"),
  names_to = "type",
  values_to = "value"
)

# Add sample group information
df_long$group <- case_when(
  grepl("INF|AS|EFF", df_long$type) ~ "WWTP",
  grepl("RU|RS|RD", df_long$type) ~ "RW",
  grepl("SW|TW", df_long$type) ~ "DW"
)
df_long$group <- factor(df_long$group, levels = c("WWTP", "RW", "DW"))
df_long$type <- factor(df_long$type, levels = c("INF", "AS", "EFF", "RU", "RS", "RD", "SW", "TW"))
df_long$ID <- factor(df_long$ID, levels = rev(c(
  "bin1180", "bin59", "bin858", "bin862", "bin896", "bin638", "bin963", "bin605", "bin447", "bin659",
  "bin689", "bin557", "bin5", "bin195", "bin1192", "bin868", "bin468", "bin140", "bin473", "bin269",
  "bin1195", "bin754", "bin68", "bin261", "bin1140", "bin149", "bin1202", "bin216", "bin1004", "bin476",
  "bin311", "bin409", "bin193", "bin1", "bin674", "bin303", "bin382", "bin349"
)))

# Define colors for bubble plot
col_bubble <- c("#c13736", "#2978b5", "#e29e4b")

# Create bubble plot
p3 <- ggplot(df_long, aes(x = type, y = ID, size = value, color = group)) +
  geom_point(
    data = subset(df_long, value != 0),
    alpha = 0.85,
    shape = 16
  ) +
  geom_point(
    data = subset(df_long, value == 0),
    alpha = 0.85,
    size = 2,
    shape = 21,
    fill = NA,
    color = "grey70"
  ) +
  scale_size(
    range = c(2, 9)
  ) +
  scale_color_manual(
    values = col_bubble,
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linetype = "dashed", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    size = guide_legend()
  )

# Read and process ARG count data for lollipop plot
data_arg <- read.csv("arg_number.csv", header = TRUE, sep = ",")
data_arg$ID <- factor(data_arg$ID, levels = rev(c(
  "bin1180", "bin59", "bin858", "bin862", "bin896", "bin638", "bin963", "bin605", "bin447", "bin659",
  "bin689", "bin557", "bin5", "bin195", "bin1192", "bin868", "bin468", "bin140", "bin473", "bin269",
  "bin1195", "bin754", "bin68", "bin261", "bin1140", "bin149", "bin1202", "bin216", "bin1004", "bin476",
  "bin311", "bin409", "bin193", "bin1", "bin674", "bin303", "bin382", "bin349"
)))

# Create ARG lollipop plot
p4 <- ggplot(data_arg, aes(ID, Value)) +
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = Value),
               linetype = "solid",
               size = 0.8,
               color = "gray60"
  ) +
  geom_point(aes(color = Value, size = Value), shape = 16, alpha = 0.95) +
  geom_text(aes(label = Value), color = "white", size = 3.5) +
  scale_size(range = c(3, 6), guide = "none") +
  scale_color_gradient(low = "#5493c4", high = "#cd5f5e", name = "ARG Count") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linetype = "dashed", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    x = "",
    y = "Number of ARGs"
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 15))

# Read and process MGE count data for lollipop plot
data_mge <- read.csv("mge_number.csv", header = TRUE, sep = ",")
data_mge$ID <- factor(data_mge$ID, levels = rev(c(
  "bin1180", "bin59", "bin858", "bin862", "bin896", "bin638", "bin963", "bin605", "bin447", "bin659",
  "bin689", "bin557", "bin5", "bin195", "bin1192", "bin868", "bin468", "bin140", "bin473", "bin269",
  "bin1195", "bin754", "bin68", "bin261", "bin1140", "bin149", "bin1202", "bin216", "bin1004", "bin476",
  "bin311", "bin409", "bin193", "bin1", "bin674", "bin303", "bin382", "bin349"
)))

# Create MGE lollipop plot
p5 <- ggplot(data_mge, aes(ID, Value)) +
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = Value),
               linetype = "solid",
               size = 0.8,
               color = "gray60"
  ) +
  geom_point(aes(color = Value, size = Value), shape = 16, alpha = 0.95) +
  geom_text(aes(label = Value), color = "white", size = 3.5) +
  scale_size(range = c(3, 8), guide = "none") +
  scale_color_gradient(low = "#e29e4b", high = "#aa3a4b") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linetype = "dashed", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    x = "",
    y = "Number of MGEs"
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100))

# Combine all plots
combined_plot <- p2 + p3 + p4 + p5 + 
  plot_layout(widths = c(1, 1.2, 0.6, 0.6)) +
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0)
  )

# Display combined plot
print(combined_plot)

# Save high-resolution figures
ggsave("ARG_MAG_Visualization.pdf", combined_plot, width = 14, height = 10, dpi = 600)
ggsave("ARG_MAG_Visualization.png", combined_plot, width = 14, height = 10, dpi = 600)
