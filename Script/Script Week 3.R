library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)

gset <- getGEO("GSE48018", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

group_info <- pData(gset)[["treatment:ch1"]]

groups <- make.names(group_info)

gset$group <- factor(groups)

nama_grup <- levels(gset$group)
print(nama_grup)

design <- model.matrix(~0 + gset$group)

colnames(design) <- levels(gset$group)

grup_vaksin <- "trivalent.influenza.vaccination"
grup_baseline <- "baseline"

contrast_formula <- paste(grup_vaksin, "-", grup_baseline)
print(paste("Kontras yang dianalisis:", contrast_formula))

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)

head(topTableResults)

probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (Vaksin Flu)",
    x = "Expression Value (log2)",
    y = "Density"
  )

umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: Baseline vs Vaccinated",
    x = "UMAP 1",
    y = "UMAP 2"
  )

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 0.1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -0.1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Vaksin Flu")

topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

topTableResults <- topTableResults %>% filter(!is.na(SYMBOL))

rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

write.csv(topTableResults, "Hasil_GSE48018_DEG.csv")