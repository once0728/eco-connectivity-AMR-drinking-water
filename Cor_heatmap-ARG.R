# --- Correlation Heatmap of ARG Types and Environmental Factors ---
# Author: once0728
# Date: 2025-10-09
# Description: Calculates Spearman correlations between ARG type abundances and environmental parameters,
# and visualizes significant relationships in a heatmap with significance labels.

# =========================================================
# 1. Package Management
# =========================================================
packages <- c("psych", "reshape2", "pheatmap", "export", "RColorBrewer")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =========================================================
# 2. Working Directory
# =========================================================
# Read data
x <- read.table(file = "ENV.csv", sep = ",", row.names = 1, header = TRUE, 
                check.names = FALSE, fileEncoding = 'GBK') 
y <- read.table(file = "ARG_Type_Abundance.csv", sep = ",", row.names = 1, header = TRUE, 
                check.names = FALSE, fileEncoding = 'GBK')

# Calculate correlation
cor <- corr.test(y, x, method = "spearman", adjust = "none") 

# Export correlation results
cmt <- cor$r
pmt <- cor$p
cmt.out <- cbind(rownames(cmt), cmt)
pmt.out <- cbind(rownames(pmt), pmt)
write.table(cmt.out, file = "cor_ARG_Type.txt", sep = "\t", row.names = FALSE)
write.table(pmt.out, file = "pvalue_ARG_Type.txt", sep = "\t", row.names = FALSE)

# Combine correlation coefficients and p-values
df <- melt(cmt, value.name = "cor")
df$pvalue <- as.vector(pmt)
head(df)
write.table(df, file = "cor-p_ARG_Type.txt", sep = "\t")

# Prepare significance labels for plotting
if (!is.null(pmt)) {
  sssmt <- pmt < 0.001
  pmt[sssmt] <- '***'
  ssmt <- pmt > 0.001 & pmt < 0.01
  pmt[ssmt] <- '**'
  smt <- pmt > 0.01 & pmt < 0.05
  pmt[smt] <- '*'
  pmt[!sssmt & !ssmt & !smt] <- ''
} else {
  pmt <- FALSE
}

# Create heatmap
p1 <- pheatmap(cmt, 
               scale = "none",
               border_color = "white",
               cluster_row = FALSE,
               cluster_col = FALSE,  
               angle_col = 90,
               fontsize = 12, 
               display_numbers = pmt, 
               fontsize_number = 10, 
               number_color = "black",
               color = colorRampPalette(colors = c("#2978b5", "white", "#c13736"))(100),
               cellwidth = 20, 
               cellheight = 20)

# Save plot
graph2pdf(file = "correlation_heatmap_ARG_Type_total.pdf", 
          font = "Times New Roman", 
          width = 10, 
          height = 8)