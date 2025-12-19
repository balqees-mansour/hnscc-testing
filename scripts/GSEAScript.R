# ================================
# Gene Set Enrichment Analysis (GSEA)
# Using ALL genes ranked by log2FoldChange
# Databases: KEGG + GO (BP, MF, CC)
# ================================

# ---- Load libraries ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)

# ---- Read DEA results ----
deg_data <- read.csv("~/Downloads/2nd trial/DE_genes_significant1.csv",
                     stringsAsFactors = FALSE)

# ---- Column names (adjust if needed) ----
gene_col  <- "hgnc_symbol"      # gene symbols
logfc_col <- "log2FoldChange"  # log2FC

# ---- Remove NA values ----
deg_data <- deg_data %>%
  filter(!is.na(.data[[logfc_col]]))

# ---- Convert SYMBOL -> ENTREZID ----
gene_map <- bitr(
  deg_data[[gene_col]],
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# ---- Merge mapping with log2FC ----
deg_mapped <- deg_data %>%
  rename(SYMBOL = hgnc_symbol) %>%  
  inner_join(gene_map, by = "SYMBOL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

# ---- Create ranked gene list ----
gene_list <- deg_mapped[[logfc_col]]
names(gene_list) <- deg_mapped$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

# ================================
# GSEA - KEGG
# ================================

gsea_kegg <- gseKEGG(
  geneList     = gene_list,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# ================================
# GSEA - GO (BP / MF / CC)
# ================================

gsea_go_bp <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_go_mf <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "MF",
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_go_cc <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "CC",
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# ================================
# Visualization: Barplots
# ================================

library(ggplot2)
library(dplyr)
library(forcats)

# 1. Data Preparation: Select top 15 pathways based on absolute NES
top_gsea <- gsea_kegg@result %>%
  arrange(desc(abs(NES))) %>%
  head(15)

# 2. Create the Bar Chart
ggplot(top_gsea, aes(
  x = fct_reorder(Description, NES), 
  y = NES,
  fill = NES
)) +
  # Use geom_col for cleaner bars
  geom_col(width = 0.75) + 
  
  # Flip coordinates for better readability of pathway names
  coord_flip() +
  
  # Set Blue for Negative and Red for Positive values
  scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    midpoint = 0
  ) + 
  
  # Clean professional theme
  theme_minimal(base_size = 14) +
  
  # Labels
  labs(
    title = "Top 15 KEGG Pathways",
    x = NULL, 
    y = "Normalized Enrichment Score (NES)",
    fill = "NES"
  ) +
  
  # Fine-tuning the appearance
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(), # Removes horizontal lines for a cleaner look
    legend.position = "right"
  )



# Prepare BP data
top_go_bp <- gsea_go_bp@result %>%
  arrange(desc(abs(NES))) %>%
  head(15)

# Plot BP
ggplot(top_go_bp, aes(x = fct_reorder(Description, NES), y = NES, fill = NES)) +
  geom_col(width = 0.7) + 
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") + 
  theme_minimal(base_size = 14) +
  labs(title = "Top 15 GO: Biological Process", x = NULL, y = "NES", fill = "NES") +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.major.y = element_blank()
  )


# Prepare MF data
# 1. Data Preparation: Filter and select top 15 MF pathways
top_go_mf <- gsea_go_mf@result %>%
  arrange(desc(abs(NES))) %>%
  head(15)

# 2. Plotting: GSEA Bar Chart for Molecular Function
ggplot(top_go_mf, aes(
  # Wrapping long pathway descriptions at 40 characters
  x = fct_reorder(str_wrap(Description, width = 40), NES), 
  y = NES, 
  fill = NES
)) +
  # Using geom_col for bar height based on NES value
  geom_col(width = 0.7) + 
  
  # Flip coordinates to make the wrapped text readable on the Y-axis
  coord_flip() +
  
  # Diverging color scale: Blue for Negative, Red for Positive, White at Zero
  scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    midpoint = 0
  ) + 
  
  # Professional clean theme
  theme_minimal(base_size = 14) +
  
  # Chart Labels
  labs(
    title = "Top 15 GO: Molecular Function", 
    x = NULL, 
    y = "Normalized Enrichment Score (NES)", 
    fill = "NES"
  ) +
  
  # Customizing text appearance
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold", color = "black", lineheight = 0.8),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank() # Removes horizontal lines for a cleaner look
  )




# Prepare CC data
top_go_cc <- gsea_go_cc@result %>%
  arrange(desc(abs(NES))) %>%
  head(15)

# Plot CC
ggplot(top_go_cc, aes(x = fct_reorder(Description, NES), y = NES, fill = NES)) +
  geom_col(width = 0.7) + 
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") + 
  theme_minimal(base_size = 14) +
  labs(title = "Top 15 GO: Cellular Component", x = NULL, y = "NES", fill = "NES") +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.major.y = element_blank()
  )


























