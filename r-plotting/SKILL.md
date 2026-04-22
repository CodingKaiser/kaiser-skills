---
name: r-single-cell-plotting
description: R plotting conventions for Seurat/scRNA-seq: UMAP square plots (aspect.ratio=1, raster=FALSE), polychrome cluster colors, DimPlot/FeaturePlot settings, saving 300 DPI PNG, volcano plots with EnhancedVolcano, statistical comparisons with rstatix/ggpubr.
---

# R Single-Cell Plotting Conventions

## Seurat / UMAP Plots

```r
cluster_colors <- unname(pals::polychrome())

# UMAP — always square, no raster
p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
             order = TRUE, repel = TRUE, cols = cluster_colors, raster = FALSE) +
  theme_bw() +
  theme(aspect.ratio = 1)
print(p)
ggsave("umap.png", p, width = 8, height = 8, dpi = 300, bg = "white")

# Feature plot
p <- FeaturePlot(obj, features = "GENE", raster = FALSE) +
  theme(aspect.ratio = 1)

# Violin (no points)
p <- VlnPlot(obj, features = "nFeature_RNA", group.by = "Sample", pt.size = 0)

# DotPlot — genes must be unique
p <- DotPlot(obj, features = unique(genes), group.by = "cell_type") + coord_flip()
```

Key rules:
- **Always** `theme(aspect.ratio = 1)` for UMAP/dim reductions
- **Always** `raster = FALSE` (prevents pixel artifacts in PDF/HTML)
- **Always** `order = TRUE, repel = TRUE` on DimPlot
- **Cluster colors**: `unname(pals::polychrome())`
- **Always** `print()` plots in Rmd chunks
- Save at **300 DPI PNG**

## Volcano Plots

```r
library(EnhancedVolcano)
EnhancedVolcano(degs,
  lab = degs$gene,
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 0.25)
```

## Statistical Comparisons (publication-ready)

```r
# Wilcoxon + Bonferroni
stat_df |>
  rstatix::wilcox_test(value ~ group) |>
  rstatix::adjust_pvalue(method = "bonferroni") |>
  rstatix::add_significance()

# Annotate ggplot directly
p + ggpubr::stat_compare_means(
  comparisons = list(c("GroupA", "GroupB")),
  method = "wilcox.test"
)
```

## General ggplot

```r
# Save any plot at 300 DPI
ggsave("plot.png", p, width = 8, height = 6, dpi = 300, bg = "white")

# Theme for publications
p + theme_bw() + theme(
  text = element_text(size = 12),
  axis.text = element_text(size = 10)
)
```
