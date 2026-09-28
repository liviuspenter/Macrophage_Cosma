library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggsignif)
library(ggrepel)

### read data

meta <- as.data.frame(readxl::read_excel("./data/20251215_metadata_combined.xlsx"))
rownames(meta) <- meta$sample
meta$condition[which(meta$condition == "vehicle_FCS")] <- "FCS"
meta$condition[which(meta$condition == "veh_AS")] <- "AS"
meta <- meta[which(meta$condition %in% c("AS", "FCS")), ]
subset <- meta$sample
samples <- paste0("./data/", subset, ".salmon")
files <- file.path(samples, "quant.sf")

names(files) <- str_replace(samples, "./data/", "") %>%
  str_replace(".salmon", "")

library(org.Hs.eg.db)
tx2gene <- read.delim("./data/tx2gene_grch38_ens94.txt")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx_id", "ensgene")], countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = T)

data <- txi$counts %>%
  round() %>%
  data.frame()

meta$patient <- factor(meta$patient, levels = unique(meta$patient))


### Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ patient + condition)
# Run analysis
dds <- DESeq(dds)

# Define contrasts for AS vs FCS
contrast_oe <- c("condition", "AS", "FCS")

# Extract results for AS vs FCS
res_tableOE <- results(dds, contrast = contrast_oe, alpha = 0.05)
padj.cutoff <- 10^-5

rld <- rlog(dds, blind = TRUE)

# Create a tibble of results
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
# Subset the tibble to keep only significant genes
sigOE <- res_tableOE_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > 1.5)

# convert normalized_counts to a data frame and transfer the row names to a new column called "gene"
normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene")

# merge together (ensembl IDs) the normalized counts data frame with a subset of the annotations in the tx2gene data frame (only the columns for ensembl gene IDs and gene symbols)
grch38annot <- tx2gene %>%
  dplyr::select(ensgene, symbol) %>%
  dplyr::distinct()

normalized_counts <- merge(normalized_counts, grch38annot, by.x = "gene", by.y = "ensgene")
normalized_counts <- normalized_counts %>%
  as_tibble()

# extract normalized expression for significant genes from the OE and control samples
norm_OEsig <- normalized_counts %>%
  filter(gene %in% sigOE$gene)
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

annotation.colors <- list(
  condition = c(AS = "#B46424", FCS = "black"),
  status = c(control = "blue", depressive = "red")
)

# Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_tableOE_tb <- res_tableOE_tb %>%
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= log2(1.5))

# Add all the gene symbols as a column from the grch38 table using bind_cols()
res_tableOE_tb <- bind_cols(res_tableOE_tb, symbol = grch38annot$symbol[match(res_tableOE_tb$gene, grch38annot$ensgene)])

# Create an empty column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% mutate(genelabels = "")

# Sort by padj values
res_tableOE_tb <- res_tableOE_tb %>% arrange(padj)

# Populate the genelabels column with contents of the gene symbols column for the first 10 rows, i.e. the top 10 most significantly expressed genes
res_tableOE_tb$genelabels[1:20] <- as.character(res_tableOE_tb$symbol[1:20])

res_tableOE_tb$significance <- "none"
res_tableOE_tb$significance[which(res_tableOE_tb$padj < 10^-5 & res_tableOE_tb$log2FoldChange > 1.5)] <- "AS"
res_tableOE_tb$significance[which(res_tableOE_tb$padj < 10^-5 & res_tableOE_tb$log2FoldChange < -1.5)] <- "FCS"

q <- ggplot() +
  ggrastr::rasterise(geom_point(data = res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj), colour = significance), size = 0.5), dpi = 600) +
  geom_hline(yintercept = -log10(10^-5), linetype = "dashed") +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +
  geom_label_repel(data = res_tableOE_tb[which(res_tableOE_tb$symbol %in% labels), ], aes(x = log2FoldChange, y = -log10(padj), label = genelabels), label.size = 0, size = 2, min.segment.length = 0) +
  scale_x_continuous("log2FC", limits = c(-10, 10)) +
  scale_y_continuous("-log10(FDR)") +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "black")) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black"),
    axis.text.y = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave(paste0("./figures/20260127_volcano_AS_FCS_combined.svg"), width = 3, height = 3, dpi = 300, plot = q)

write.table(res_tableOE_tb[which(res_tableOE_tb$significance != "none"), c("symbol", "log2FoldChange", "padj")],
  file = "./data/20260729_diff_genes.csv", sep = "\t", dec = ".", row.names = T
)


### visualize in heatmap

colnames(norm_OEsig)[2:17] <- rownames(meta)

col_fun <- circlize::colorRamp2(breaks = seq(-3, 3, 6 / 8), colors = BuenColors::jdb_palette(name = "brewer_yes"))
ha <- columnAnnotation(
  condition = meta$condition, col = list("condition" = c("AS" = "#B46424", "FCS" = "lightgrey")),
  simple_anno_size = unit(5, "pt"), border = T
)

labels <- c(
  "PLIN2", "FABP4", "PDK4", "DHRS9", "SDS", "PDK4", "DHRS9", "RRM2", "FN1", "SPON2", "ADAMTS15", "MMP12",
  "SCIN", "C1QB", "A2M", "CXCL16", "PADI2", "PLK1", "CDK1", "TPX2", "RRM2", "MKI67",
  "CCL2", "CXCL8", "IL1B", "CD300E", "PLXNA2", "NQO1", "PDPN", "S100A8", "SMAD6", "CYP1B1", "PID1", "SGPP2"
)

norm_OEsig <- as.data.frame(norm_OEsig)
rownames(norm_OEsig) <- norm_OEsig$symbol


ha_row <- rowAnnotation(foo = anno_mark(
  at = match(labels, norm_OEsig$symbol),
  side = "left", labels = labels, labels_gp = gpar(fontsize = 8)
))
svglite::svglite(paste0("./figures/20260127_heatmap_AS_FCS_combined.svg"), width = 3.5, height = 5)
Heatmap(t(scale(t(norm_OEsig[2:17]))),
  top_annotation = ha, show_column_dend = F, show_row_dend = F,
  show_column_names = F, border = T, left_annotation = ha_row,
  column_split = factor(meta$condition, levels = c("AS", "FCS")), col = col_fun
)
dev.off()


### correlation between conditions
library(reshape2)

expr.mat <- normalized_counts[, c(make.names(meta$sample[which(meta$condition == "AS")]), make.names(meta$sample[which(meta$condition == "FCS")]))]
rownames(expr.mat) <- normalized_counts$gene

cormat <- cor(expr.mat[sigOE$gene, ])

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)
# Heatmap

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "firebrick", mid = "white",
    midpoint = 0.55, limit = c(0, 1), space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(
    angle = 45, vjust = 1,
    size = 12, hjust = 1
  )) +
  coord_fixed()

ggheatmap +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  ))
ggsave("./figures/20260927_correlation_between_conditions.svg", width = 4.5, height = 3.5)

### visualize genes of interest

normalized_counts.collapsed <- normalized_counts %>%
  group_by(symbol) %>%
  summarize(across(where(is.numeric), sum, na.rm = T))

expr.mat <- t(normalized_counts.collapsed) %>% as.data.frame()
colnames(expr.mat) <- expr.mat[1, ]
expr.mat <- as.data.frame(sapply(expr.mat[-1, ], as.numeric))
expr.mat$condition <- meta$condition
expr.mat$patient <- meta$patient

ggplot(expr.mat, aes(x = CD300A, y = CD300E, color = condition)) +
  geom_point() +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260127_CD300A_CD300E_combined.svg", width = 2.2, height = 2.2)

ggplot(expr.mat, aes(x = C1QB, y = IL1B, color = condition)) +
  geom_point() +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260127_C1QB_IL1B_combined.svg", width = 2.2, height = 2.2)

ggplot(expr.mat, aes(x = PDK4, y = NQO1, color = condition)) +
  geom_point() +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260127_PDK4_NQO1_combined.svg", width = 2.2, height = 2.2)

ggplot(expr.mat, aes(x = FN1, y = PDPN, color = condition)) +
  geom_point() +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260127_FN1_PDPN_combined.svg", width = 2.2, height = 2.2)

ggplot(expr.mat, aes(x = condition, y = CD300A)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_CD300A.svg", width = 1.3, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = C1QB)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_C1QB.svg", width = 1.3, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = IL1B)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_IL1B.svg", width = 1.3, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = PDK4)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_PDK4.svg", width = 1.3, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = NQO1)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_NQO1.svg", width = 1.3, height = 2.5)


ggplot(expr.mat, aes(x = condition, y = FN1)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_FN1.svg", width = 1.3, height = 2.5)


ggplot(expr.mat, aes(x = condition, y = PDPN)) +
  geom_signif(comparisons = list(c("AS", "FCS")), textsize = 3) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, color = "black") +
  geom_line(aes(group = patient)) +
  geom_point(aes(color = condition), size = 0.5) +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figures/20260202_FN1.svg", width = 1.3, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = TNF, color = condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  geom_signif(comparisons = list(c("AS", "FCS")), color = "black") +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./figures/20260127_TNF_combined.svg", width = 1.5, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = IL1B, color = condition)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  geom_signif(comparisons = list(c("AS", "FCS")), color = "black") +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./figures/20260127_IL1B_combined.svg", width = 1.5, height = 2.5)

ggplot(expr.mat, aes(x = condition, y = IL6, color = condition)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  geom_signif(comparisons = list(c("AS", "FCS")), color = "black") +
  scale_color_manual(values = c("none" = "lightgrey", AS = "#B46424", FCS = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./figures/20260127_IL6_combined.svg", width = 1.5, height = 2.5)
