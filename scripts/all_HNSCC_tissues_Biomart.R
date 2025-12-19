############### 1. Load Required Libraries ###############
# Data handling and visualization

library(ggplot2)
library(dplyr)

# TCGA data handling
library(curatedTCGAData)
library(TCGAbiolinks)
library(SummarizedExperiment)

# Differential Expression
library(DESeq2)

library(sva)

# Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

############### 2. Data Retrieval ###############

clinical_hnsc <- GDCquery_clinic("TCGA-HNSC")

# Build query for TCGA-COAD RNA-Seq data
query_TCGA <- GDCquery(
  project = 'TCGA-HNSC',
  data.category = c('Transcriptome Profiling'),
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  access = 'open'
)


# Download the data
GDCdownload(query_TCGA, directory = "/home/balqees/Desktop/HNSCC/TCGA_updated analysis")

# Prepare the data
hnsc.tcga.data <- GDCprepare(query_TCGA, 
                             summarizedExperiment = TRUE,
                             directory = "/home/balqees/Desktop/HNSCC/TCGA_updated analysis")

#hnsc_normal <- hnsc.tcga.data[, hnsc.tcga.data$definition == "Solid Tissue Normal"]



table(hnsc.tcga.data@colData$definition)
tissues <- as.data.frame(table(hnsc.tcga.data@colData@listData[["tissue_or_organ_of_origin"]]))

clinical <- as.data.frame(colData(hnsc.tcga.data))
cases <- read.csv("/home/balqees/Desktop/HNSCC/TCGA_updated analysis/cases.tsv", sep = '\t')


final_barcodes <- clinical[which(cases$submitter_id %in%  clinical$patient),]

table(final_barcodes$definition)
barcodes <- final_barcodes$barcode


# Build query for TCGA-COAD RNA-Seq data
query_TCGA_subset <- GDCquery(
  project = 'TCGA-HNSC',
  data.category = c('Transcriptome Profiling'),
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  barcode = barcodes , 
  access = 'open'
)

# Download the data
#GDCdownload(query_TCGA_subset, directory = "/home/balqees/Desktop/HNSCC/TCGA_updated analysis")

# Prepare the data
hnsc.tcga.data <- GDCprepare(query_TCGA_subset, 
                             summarizedExperiment = TRUE,
                             directory = "/home/balqees/Desktop/HNSCC/TCGA_updated analysis")



table(hnsc.tcga.data$definition)

# Filter out metastatic samples
hnsc_no_metastatic <- hnsc.tcga.data[, hnsc.tcga.data$definition != "Metastatic"]

# Check the definition table after removal
table(hnsc_no_metastatic$definition)

# Extract count data
counts_data <- assay(hnsc_no_metastatic)

#write.csv(counts_data , "counts_tcga_data.csv")

# # Create DESeq2 object
# dds <- DESeqDataSetFromMatrix(countData = counts_data,
#                               colData = colData(hnsc_no_metastatic),
#                               design = ~ 1) # Design can be adjusted later for differential analysis
# 
# 
# 
# # Perform VST normalization for PCA
# vsd <- vst(dds, blind = TRUE)
# 
# # Perform PCA
# pca_data <- plotPCA(vsd, intgroup = "definition", returnData = TRUE)
# percentVar <- round(100 * attr(pca_data, "percentVar"))
# 
# # Plot PCA
# ggplot(pca_data, aes(PC1, PC2, color=definition)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed() +
#   ggtitle("PCA of TCGA-HNSC RNA-Seq Data")
# 
# # Save PCA plot (optional)
# ggsave("pca_plot.png", width = 8, height = 6, units = "in")



hnsc_no_metastatic_counts <- assay(hnsc_no_metastatic)

hnsc_no_metastatic_counts[,"TCGA-CV-7446-01A-11R-2232-07"]

pheno.data <- colData(hnsc_no_metastatic)
pheno.data <- as.data.frame(pheno.data)
pheno.data <- pheno.data[,c(19,24)]
#write.csv(pheno.data, "sample_sheet2.csv")

############### 4. DESeq2 Analysis ###############
# Create DESeq dataset
pheno.data$sample_type = factor(pheno.data$tissue_type, 
                                levels = c("Normal", "Tumor"))

dds<- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData =  colData(hnsc_no_metastatic),
                              design = ~ sample_type)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]



############### 5. Batch Effect Correction ###############
# Perform SVA analysis
# Normalize and filter
dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)

idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ sample_type, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

# Create data frame for ggplot2
plot_data <- data.frame(
  sv1 = svseq$sv[, 1],
  sv2 = svseq$sv[, 2],
  tissue = dds$tissue_or_organ_of_origin,
  stringsAsFactors = FALSE
)

# Manual pivot to long format (avoiding tidyr dependency)
plot_data_long <- data.frame(
  tissue = rep(plot_data$tissue, 2),
  value = c(plot_data$sv1, plot_data$sv2),
  sv_component = rep(c("sv1", "sv2"), each = nrow(plot_data)),
  stringsAsFactors = FALSE
)

# Enhanced boxplot with better styling
ggplot(plot_data_long, aes(x = tissue, y = value, fill = sv_component)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.6) +
  facet_wrap(~sv_component, scales = "free_y", ncol = 1, 
             labeller = labeller(sv_component = c("sv1" = "SV1 by tissue", 
                                                  "sv2" = "SV2 by tissue"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("sv1" = "#56B4E9", "sv2" = "#E69F00")) +
  labs(x = "Tissue Type", y = "SV Component Value")


# Create a data frame for SVA results
sva_df <- data.frame(
  SV1 = svseq$sv[,1],
  SV2 = svseq$sv[,2],
  condition = colData(dds)$sample_type
)
# Create scatter plot of SV1 vs SV2
ggplot(sva_df, aes(x = SV1, y = SV2, color = condition)) +
  geom_point(size = 6, alpha = 0.8) +
  scale_color_manual(values = c("Primary Tumor" = "#FF9999", 
                                "Solid Tissue Normal" = "#66CCCC")) +
  xlab("Surrogate Variable 1") +
  ylab("Surrogate Variable 2") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    legend.position = "right",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Save the plot
ggsave("sva_plotl.png", width = 10, height = 6, dpi = 300)


# Add surrogate variables to colData
colData(dds)$SV1 <- svseq$sv[,1]
colData(dds)$SV2 <- svseq$sv[,2]

# Update design formula
design(dds) <- ~ SV1 + SV2 + sample_type

# Run DESeq
dds <- DESeq(dds)

plotDispEsts(dds, main="Dispersion estimates and adjustments")

plotMA(dds, main="MA plot: FDR<0.1 for any fold change")

# Calculate results for 2-fold change and FDR<0.05 
# Note that in log2-scale the 2-fold change corresponds to the lfcThreshold=1
FDR_0.05_2x_FC <- results(dds, lfcThreshold=1.5, alpha=0.05)

# Explore results
summary(FDR_0.05_2x_FC)

plotMA(dds, lfcThreshold=1, alpha=0.05, main="MA plot: FDR<0.05 for at least 2x change")

# Extract results to data frame
result.df <- as.data.frame(FDR_0.05_2x_FC)

# Input
ens_ids <- rownames(result.df)

# 1) Strip version suffix
ens_stable <- sub("\\.\\d+$", "", ens_ids)  # ENSG... without version [16]

# 2) Query Ensembl via biomaRt
suppressPackageStartupMessages({
  library(biomaRt)
})

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  # connects to current Ensembl [3]
anno <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
  filters = "ensembl_gene_id",
  values = ens_stable,
  mart = mart
)

# 3) Join back to original order and include versioned IDs
res_sym <- data.frame(ensembl_with_version = ens_ids, ensembl_gene_id = ens_stable)
res_sym <- merge(res_sym, anno, by = "ensembl_gene_id", all.x = TRUE)
res_sym <- res_sym[match(ens_stable, res_sym$ensembl_gene_id), c("ensembl_with_version","ensembl_gene_id","hgnc_symbol","gene_biotype","description")]
colnames(res_sym)[1] <- "symbol"

res_sym
result.df$ensembl_gene_id <- rownames(result.df)

dd <- merge(res_sym, result.df, by = "ensembl_gene_id", all.x = TRUE)

result <- res_sym %>%
  left_join(result.df, by = "ensembl_gene_id")

res_sym <- (na.omit(res_sym))
# 

library(dplyr)
res_sym <- res_sym %>% mutate(ensembl_gene_id = sub("\\.\\d+$", "", ensembl_gene_id))
result.df <- result.df %>% mutate(ensembl_gene_id = sub("\\.\\d+$", "", ensembl_gene_id))
intersect_n <- dplyr::n_distinct(intersect(res_sym$ensembl_gene_id, result.df$ensembl_gene_id))
result <- res_sym %>% left_join(result.df, by = "ensembl_gene_id")

result <- na.omit(result)

result$hgnc_symbol
result <- result[,-1]
# Select columns to save and order by p-value (and fold change)
result <- result %>% 
  dplyr::select(ensembl_gene_id, hgnc_symbol, baseMean, log2FoldChange,  padj) %>% 
  arrange(padj, desc(log2FoldChange))
dim(result)



mu_gene_updated <- read.csv("/home/balqees/Desktop/HNSCC/TCGA_updated analysis/MU gene list-string.csv")
new_list_results <- as.data.frame(result[result$hgnc_symbol %in% mu_gene_updated$gene_id ,])

new_filter_all <- new_list_results[which(new_list_results$padj <= 0.05),]
write.csv(new_filter_all, "filtered_mitgenes_from_TCGAdata.csv")


# i.e. the thresholds that could be adjusted for specific study purpose
degs.df <- result %>% 
  filter(padj <= 10e-5, abs(log2FoldChange) >= 5)

# Save Ensemble gene IDs
library(stringr)
degs_ens_id.df = degs.df %>% 
  mutate(ens_gene_id = str_remove(gene_id, '\\..*')) %>% 
  dplyr::select(ens_gene_id)
degs_ens_file <- file.path(output_folder,"top_DEGs_ens_ids.txt")
write.table(degs_ens_id.df, file=degs_ens_file, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

# Save genes names (NB: mix of Symbols and Ensemble ids!)
degs_gene_name.df = degs.df %>% 
  dplyr::select(hgnc_symbol)
degs_names_file <- file.path(output_folder,"top_DEGs_names.txt")
write.table(degs_gene_name.df, file=degs_names_file, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

# Clean-up
rm(excluded_genes, results_file, output_folder, gene_names.df, 
   degs.df, degs_ens_id.df, degs_ens_file, degs_gene_name.df, degs_names_file)



library(EnhancedVolcano)
#?EnhancedVolcano

EnhancedVolcano(result,
                lab = result$hgnc_symbol,
                x = "log2FoldChange",
                y = "padj", 
                title="DESeq2 results",
                subtitle="Enhanced Volcano Plot")















#==============================================================================
#==============================================================================
#-----------01# BoxPLot for the 20 upregulated genes in tumor vs normal -------
#------------------------------------------------------------------------------
#==============================================================================
#==============================================================================

library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'
# Order by log2FoldChange
res_upregulated <- res_all_tiisues[order(res_all_tiisues$padj, decreasing = FALSE),]
res_upregulated <- res_upregulated[which(res_upregulated$padj < 0.05 & res_upregulated$Significance == "Upregulated") ,]

# Select top 18 upregulated genes
top_upregulated <- head(res_upregulated, 25)

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(top_upregulated),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(pheno.data, by = c("sample" = "ids"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(top_upregulated) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = sample_type, y = log2(count + 1), fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("Normal" = "lightblue", "Tumor" = "indianred")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = sample_type[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_upregulated_genes_HNSCC.png", width = 18, height = 14, dpi = 300)

write.csv(top_upregulated, "top_20_upregulated_genes_HNSCC.csv")
#==============================================================================
#----- #02 BoxPLot for the top 20  downregulated genes in tumor vs normal -----
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'


# Order by log2FoldChange
res_downregulated <- res_all_tiisues[order(res_all_tiisues$padj, decreasing = FALSE),]
res_downregulated <- res_downregulated[which(res_downregulated$padj < 0.05 & res_downregulated$Significance == "Downregulated") ,]

# Select top 18 upregulated genes
top_downregulated <- head(res_downregulated, 25)

# Extract normalized counts for these genes
normalized_counts_d <- counts(dds, normalized=TRUE)[rownames(top_downregulated),]
pheno.data$ids <- rownames(pheno.data)
# Prepare data for plotting
plot_data <- normalized_counts_d %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(pheno.data, by = c("sample" = "ids"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(top_downregulated) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = sample_type, y = log2(count + 1), fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("Normal" = "lightblue", "Tumor" = "indianred")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = sample_type[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_downregulated_genes.png", width = 18, height = 14, dpi = 300)
write.csv(top_downregulated,"Downregulated_20genes_HNSCC.csv")
#===========================================================
## GSEA
#===========================================================
library(clusterProfiler)
library(org.Hs.eg.db)
# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(top_upregulated)
gene_list <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filter out any NA values that may result from conversion
gene_list <- gene_list[!is.na(gene_list$ENTREZID), "ENTREZID"]

# Perform KEGG enrichment analysis
ekegg <- enrichKEGG(gene         = gene_list,
                    organism     = 'hsa', # Human
                    pvalueCutoff = 0.05)

# Convert enrichment results to a data frame
ekegg_df <- as.data.frame(ekegg@result)

# Sort the data frame by raw p-value in ascending order
ekegg_df <- ekegg_df[order(ekegg_df$pvalue), ]

# Select the top 10 significant pathways based on p-value
ekegg_df <- ekegg_df[1:10, ]

# Calculate Rich Factor
ekegg_df$RichFactor <- ekegg_df$Count / as.numeric(sapply(strsplit(ekegg_df$GeneRatio, "/"), `[`, 2))

# Create a lollipop plot using ggplot2
library(ggplot2)
library(forcats)

ggplot(ekegg_df, aes(x = RichFactor, y = fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = pvalue, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("KEGG Enrichment Lollipop PlotTop 20 Upregulated Genes)") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Increase title size
    axis.title.x = element_text(size = 14),                  # Increase x-axis title size
    axis.text.x = element_text(size = 14),                   # Increase x-axis text size
    axis.text.y = element_text(size = 14),                   # Increase y-axis text size
    legend.title = element_text(size = 14),                  # Increase legend title size
    legend.text = element_text(size = 12)                    # Increase legend text size
  )
# Save the plot
ggsave("lollipop_KEGG10_uppathways_pvalue.png", width = 15, height = 10, dpi = 300)


# down 20 
library(clusterProfiler)
library(org.Hs.eg.db)
# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(top_downregulated)
gene_list <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filter out any NA values that may result from conversion
gene_list <- gene_list[!is.na(gene_list$ENTREZID), "ENTREZID"]

# Perform KEGG enrichment analysis
ekegg <- enrichKEGG(gene         = gene_list,
                    organism     = 'hsa', # Human
                    pvalueCutoff = 0.05)

# Convert enrichment results to a data frame
ekegg_df <- as.data.frame(ekegg@result)

# Sort the data frame by raw p-value in ascending order
ekegg_df <- ekegg_df[order(ekegg_df$pvalue), ]

# Select the top 10 significant pathways based on p-value
ekegg_df <- ekegg_df[1:10, ]

# Calculate Rich Factor
ekegg_df$RichFactor <- ekegg_df$Count / as.numeric(sapply(strsplit(ekegg_df$GeneRatio, "/"), `[`, 2))

# Create a lollipop plot using ggplot2
library(ggplot2)
library(forcats)
ggplot(ekegg_df, aes(x = RichFactor, y = fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = pvalue, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("KEGG Enrichment Lollipop Plot (Top 20 Downregulated Genes)") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),      # Increase title size
    axis.title.x = element_text(size = 16),                  # Increase x-axis title size
    axis.text.x = element_text(size = 14),                   # Increase x-axis text size
    axis.text.y = element_text(size = 14),                   # Increase y-axis text size
    legend.title = element_text(size = 14),                  # Increase legend title size
    legend.text = element_text(size = 12)                    # Increase legend text size
  )


# Save the plot
ggsave("lollipop_KEGG10down_pathways_pvalue.png", width = 15, height = 10, dpi = 300)







