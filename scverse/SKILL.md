---
name: scverse
description: Single-cell analysis patterns with AnnData and Scanpy. Use when working with single-cell RNA-seq data, creating or manipulating AnnData objects, performing quality control, normalization, dimensionality reduction, clustering, or visualization of single-cell datasets.
---

# scverse: AnnData and Scanpy Guide

Patterns for single-cell analysis using the scverse ecosystem, focusing on the AnnData data structure and Scanpy analysis workflows.

## Quick Reference

```python
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

# Read/write data
adata = sc.read_h5ad("data.h5ad")
adata.write_h5ad("output.h5ad")

# Basic QC and preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Standard workflow
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Plotting
sc.pl.umap(adata, color=["leiden", "gene_name"])
```

## AnnData Structure

AnnData is the core data structure for single-cell data:

```
adata
├── .X                  # Expression matrix (cells × genes)
├── .obs                # Cell metadata (DataFrame)
├── .var                # Gene metadata (DataFrame)
├── .uns                # Unstructured data (dict)
├── .obsm               # Cell embeddings (e.g., PCA, UMAP)
├── .varm               # Gene embeddings
├── .obsp               # Cell-cell graphs (e.g., neighbors)
├── .varp               # Gene-gene relationships
└── .layers             # Alternative expression matrices
```

### Creating AnnData Objects

```python
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# From dense matrix
X = np.random.poisson(1, size=(100, 2000))  # 100 cells, 2000 genes
adata = ad.AnnData(X)

# With metadata
adata = ad.AnnData(
    X=X,
    obs=pd.DataFrame({"cell_type": ["A"] * 50 + ["B"] * 50}, index=[f"cell_{i}" for i in range(100)]),
    var=pd.DataFrame({"gene_name": [f"gene_{i}" for i in range(2000)]}, index=[f"gene_{i}" for i in range(2000)])
)

# From sparse matrix (recommended for large datasets)
X_sparse = csr_matrix(X)
adata = ad.AnnData(X_sparse)

# From DataFrame
df = pd.DataFrame(X, index=[f"cell_{i}" for i in range(100)], columns=[f"gene_{i}" for i in range(2000)])
adata = ad.AnnData(df)
```

### Accessing and Modifying Data

```python
# Expression matrix
adata.X                    # Full matrix
adata.X[0, :]              # First cell
adata[:, "gene_1"].X       # Single gene across all cells

# Metadata
adata.obs["cluster"]       # Cell annotations
adata.var["highly_variable"]  # Gene annotations
adata.obs_names            # Cell IDs
adata.var_names            # Gene names

# Shape
adata.n_obs                # Number of cells
adata.n_vars               # Number of genes
adata.shape                # (n_cells, n_genes)

# Add new annotations
adata.obs["new_column"] = values
adata.var["is_mito"] = adata.var_names.str.startswith("MT-")
```

### Subsetting

```python
# By index
adata_subset = adata[0:100, :]           # First 100 cells
adata_subset = adata[:, 0:500]           # First 500 genes

# By boolean mask
mask = adata.obs["cell_type"] == "T cell"
adata_tcells = adata[mask, :]

# By names
cells_of_interest = ["cell_1", "cell_2", "cell_3"]
adata_subset = adata[cells_of_interest, :]

genes_of_interest = ["CD3D", "CD4", "CD8A"]
adata_subset = adata[:, genes_of_interest]

# Combined
adata_subset = adata[mask, genes_of_interest]

# Copy vs view (important!)
adata_view = adata[mask, :]              # View (linked to original)
adata_copy = adata[mask, :].copy()       # Independent copy
```

### Layers

```python
# Store raw counts before normalization
adata.layers["counts"] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Access layers
adata.layers["counts"]         # Raw counts
adata.X                        # Normalized data

# Use layers in plotting
sc.pl.violin(adata, "gene_name", layer="counts")
```

### Concatenation

```python
# Concatenate samples
adata_combined = ad.concat(
    [adata1, adata2, adata3],
    join="outer",              # Keep all genes
    label="sample",            # Add sample label
    keys=["sample1", "sample2", "sample3"]
)

# Merge specific batches
adata_combined = ad.concat(
    {"batch1": adata1, "batch2": adata2},
    label="batch"
)
```

## I/O Patterns

### Reading Data

```python
# H5AD (native format)
adata = sc.read_h5ad("data.h5ad")

# 10x Genomics
adata = sc.read_10x_mtx("filtered_feature_bc_matrix/")
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")

# CSV/TSV
adata = sc.read_csv("expression.csv")
adata = sc.read_text("expression.tsv", delimiter="\t")

# Loom
adata = sc.read_loom("data.loom")
```

### Writing Data

```python
# H5AD (recommended)
adata.write_h5ad("output.h5ad")
adata.write_h5ad("output.h5ad", compression="gzip")

# CSV (dense, for small datasets)
adata.to_df().to_csv("expression.csv")

# Loom
adata.write_loom("output.loom")
```

## Preprocessing Workflow

### Quality Control

```python
# Calculate QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")  # Mitochondrial genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))  # Ribosomal

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Adds to adata.obs:
# - n_genes_by_counts: genes detected per cell
# - total_counts: total UMIs per cell
# - pct_counts_mt: % mitochondrial reads
# - pct_counts_ribo: % ribosomal reads

# Visualize QC metrics
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt")
```

### Filtering

```python
# Filter cells
sc.pp.filter_cells(adata, min_genes=200)     # At least 200 genes
sc.pp.filter_cells(adata, max_genes=5000)    # Remove doublets

# Filter by QC metrics
adata = adata[adata.obs["pct_counts_mt"] < 20, :]  # Max 20% mito

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)       # Gene in at least 3 cells

# Combined filtering
adata = adata[
    (adata.obs["n_genes_by_counts"] > 200) &
    (adata.obs["n_genes_by_counts"] < 5000) &
    (adata.obs["pct_counts_mt"] < 20),
    :
].copy()
```

### Normalization

```python
# Store raw counts first
adata.layers["counts"] = adata.X.copy()

# Or use raw slot
adata.raw = adata

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
sc.pp.log1p(adata)

# Alternative: scran-style normalization (for batch effects)
# Requires the scran package
```

### Highly Variable Genes

```python
# Identify HVGs
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,          # Select top 2000 HVGs
    flavor="seurat_v3",        # Method (seurat, seurat_v3, cell_ranger)
    batch_key="sample"         # Optional: per-batch HVG selection
)

# Adds to adata.var:
# - highly_variable: boolean
# - means, dispersions, dispersions_norm

# Visualize
sc.pl.highly_variable_genes(adata)

# Subset to HVGs (optional)
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
```

### Scaling

```python
# Scale to unit variance (optional, mainly for PCA)
sc.pp.scale(adata, max_value=10)  # Clip to avoid outlier influence

# Note: scaling modifies X, consider using layers
adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X
```

## Dimensionality Reduction

### PCA

```python
# Run PCA on HVGs
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

# Results stored in:
# - adata.obsm["X_pca"]: cell embeddings
# - adata.varm["PCs"]: gene loadings
# - adata.uns["pca"]["variance_ratio"]: explained variance

# Visualize
sc.pl.pca(adata, color="cell_type")
sc.pl.pca_variance_ratio(adata, n_pcs=50)

# Determine number of PCs to use
sc.pl.pca_variance_ratio(adata, log=True)
```

### Neighborhood Graph

```python
# Build k-NN graph (required for UMAP, clustering)
sc.pp.neighbors(
    adata,
    n_neighbors=15,      # Number of neighbors
    n_pcs=30,            # PCs to use
    metric="euclidean"
)

# Results stored in:
# - adata.obsp["distances"]: distance matrix
# - adata.obsp["connectivities"]: adjacency matrix
# - adata.uns["neighbors"]: parameters
```

### UMAP

```python
# Compute UMAP
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# Results in adata.obsm["X_umap"]

# Visualize
sc.pl.umap(adata, color=["cell_type", "n_genes_by_counts"])

# Save high-resolution figure
sc.pl.umap(adata, color="cell_type", save="_celltype.pdf")
```

### t-SNE

```python
# Compute t-SNE (slower than UMAP)
sc.tl.tsne(adata, n_pcs=30, perplexity=30)

# Results in adata.obsm["X_tsne"]

sc.pl.tsne(adata, color="cell_type")
```

## Clustering

### Leiden Clustering

```python
# Leiden clustering (recommended over Louvain)
sc.tl.leiden(
    adata,
    resolution=0.5,      # Lower = fewer clusters
    key_added="leiden"   # Column name in adata.obs
)

# Try multiple resolutions
for res in [0.3, 0.5, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")

# Visualize
sc.pl.umap(adata, color=["leiden_0.3", "leiden_0.5", "leiden_0.8"])
```

### Louvain Clustering

```python
# Louvain (alternative to Leiden)
sc.tl.louvain(adata, resolution=0.5)

sc.pl.umap(adata, color="louvain")
```

## Differential Expression

### Rank Genes by Groups

```python
# Find marker genes per cluster
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",     # wilcoxon, t-test, logreg
    pts=True               # Include fraction of cells expressing
)

# Results in adata.uns["rank_genes_groups"]

# Visualize
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5)

# Extract as DataFrame
markers = sc.get.rank_genes_groups_df(adata, group="0")  # Cluster 0 markers
markers_all = sc.get.rank_genes_groups_df(adata, group=None)  # All clusters
```

### Filter Markers

```python
# Get significant markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
significant = markers[
    (markers["pvals_adj"] < 0.05) &
    (markers["logfoldchanges"] > 1)
]
```

## Visualization

### UMAP/t-SNE Plots

```python
# Basic UMAP
sc.pl.umap(adata, color="leiden")

# Multiple annotations
sc.pl.umap(adata, color=["leiden", "cell_type", "sample"])

# Gene expression
sc.pl.umap(adata, color=["CD3D", "CD4", "CD8A"], cmap="viridis")

# Customize appearance
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    legend_fontsize=8,
    frameon=False,
    title="Cell Clusters"
)
```

### Dot Plots

```python
# Marker gene expression by cluster
marker_genes = ["CD3D", "CD4", "CD8A", "MS4A1", "CD14", "FCGR3A"]
sc.pl.dotplot(adata, marker_genes, groupby="leiden")

# With dendrogram
sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="leiden",
    dendrogram=True,
    standard_scale="var"  # Scale per gene
)
```

### Heatmaps

```python
# Expression heatmap
sc.pl.heatmap(
    adata,
    marker_genes,
    groupby="leiden",
    cmap="viridis",
    dendrogram=True
)

# Rank genes heatmap
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby="leiden")
```

### Violin Plots

```python
# Gene expression by cluster
sc.pl.violin(adata, ["CD3D", "CD4"], groupby="leiden")

# Stacked violin
sc.pl.stacked_violin(adata, marker_genes, groupby="leiden")
```

### Matrix Plots

```python
# Mean expression per group
sc.pl.matrixplot(
    adata,
    marker_genes,
    groupby="leiden",
    dendrogram=True,
    cmap="Blues"
)
```

## Best Practices

### Workflow Checklist

- [ ] Store raw counts in `adata.layers["counts"]` or `adata.raw` before normalization
- [ ] Check QC metrics and filter appropriately for your tissue/technology
- [ ] Use `flavor="seurat_v3"` for HVG selection with raw counts
- [ ] Choose appropriate number of PCs (elbow in variance ratio plot)
- [ ] Try multiple clustering resolutions and validate biologically
- [ ] Use Leiden over Louvain (better performance)
- [ ] Save intermediate results with `adata.write_h5ad()`

### Memory Management

```python
# Use sparse matrices
from scipy.sparse import csr_matrix
if not isinstance(adata.X, csr_matrix):
    adata.X = csr_matrix(adata.X)

# Process in chunks for very large datasets
sc.settings.n_jobs = 4  # Parallel processing

# Clear cache
import gc
gc.collect()
```

### Reproducibility

```python
# Set random seed
import numpy as np
np.random.seed(42)

# For UMAP
sc.tl.umap(adata, random_state=42)

# For clustering
sc.tl.leiden(adata, random_state=42)
```

## Resources

- [AnnData Documentation](https://anndata.readthedocs.io/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [scverse Ecosystem](https://scverse.org/)
- [Best Practices for Single-Cell Analysis](https://www.sc-best-practices.org/)

## Common Issues

### Memory errors with large datasets

```python
# Use backed mode for very large files
adata = sc.read_h5ad("large_data.h5ad", backed="r")

# Or use sparse matrices throughout
adata.X = csr_matrix(adata.X)
```

### Gene names not matching

```python
# Check if using symbols vs Ensembl IDs
print(adata.var_names[:5])

# Convert if needed (requires biomart or similar)
# Or check var columns for alternative names
print(adata.var.columns)
```

### Plotting issues

```python
# Set figure parameters
sc.settings.set_figure_params(dpi=100, facecolor="white")

# Save figures to specific directory
sc.settings.figdir = "./figures/"
```
