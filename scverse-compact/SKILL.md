---
name: scverse-compact
description: Compact single-cell analysis patterns with AnnData and Scanpy. Essential patterns for QC, preprocessing, clustering, and visualization.
---

# scverse: AnnData & Scanpy Essentials

## AnnData Structure

```
adata.X          # Expression matrix (cells × genes)
adata.obs        # Cell metadata (DataFrame)
adata.var        # Gene metadata (DataFrame)
adata.uns        # Unstructured data (dict)
adata.obsm       # Embeddings (PCA, UMAP)
adata.layers     # Alternative matrices (raw counts)
```

## Standard Workflow
```python
import scanpy as sc

# Read data
adata = sc.read_h5ad("data.h5ad")

# QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()

# Normalize (store raw first)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG, PCA, neighbors, UMAP, clustering
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Markers
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
```

## Key Patterns

**Subsetting**: `adata[mask, genes].copy()` - always `.copy()` for independent object

**Concatenation**: `ad.concat([adata1, adata2], label="sample", keys=["s1", "s2"])`

**Plotting**:
- `sc.pl.umap(adata, color=["leiden", "gene"])`
- `sc.pl.dotplot(adata, marker_genes, groupby="leiden")`
- `sc.pl.violin(adata, ["gene1", "gene2"], groupby="leiden")`

## Critical Rules
- Store raw counts before normalization: `adata.layers["counts"] = adata.X.copy()`
- Use `.copy()` after subsetting to avoid views
- Leiden > Louvain for clustering
- Set `random_state=42` for reproducibility
