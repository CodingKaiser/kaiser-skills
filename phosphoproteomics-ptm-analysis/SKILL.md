# Phosphoproteomics PTM Analysis Skill

This skill provides guidance for analyzing phosphoproteomics data using R, focusing on differential PTM abundance (DPA), differential PTM usage (DPU), and kinase enrichment analysis workflows.

## Overview

Phosphoproteomics analysis examines post-translational modifications (PTMs), specifically phosphorylation, to understand cell signaling changes. This skill covers:

- **DPA (Differential PTM Abundance)**: Raw phosphosite signal changes
- **DPU (Differential PTM Usage)**: Protein-normalized stoichiometry changes
- **CorrectFirst**: Alternative protein correction before statistical modeling
- **KSEA**: Kinase-Substrate Enrichment Analysis (gene + position based)
- **PTM-SEA**: PTM Signature Enrichment Analysis (flanking sequence based, species-invariant)

## R Package Ecosystem

### Core Proteomics Analysis

```r
# Statistical framework for differential expression
library(prolfqua)       # LFQData, Contrasts, statistical modeling
library(prolfquapp)     # High-level DEA pipeline, configuration management
library(prophosqua)     # Phosphoproteomics-specific (DPA/DPU calculations)

# Data readers (search engine specific)
library(prolfquappPTMreaders)  # Spectronaut BGS, FragPipe, MaxQuant support
```

**Key prolfqua patterns:**
```r
# R6 object creation
lfq_data <- prolfqua::LFQData$new(lfq_config, data)

# Model building
strategy <- prolfqua::strategy_lm()
model_result <- prolfqua::build_model(lfq_data, strategy)

contrasts <- prolfqua::Contrasts$new(model_result, contrast_def)

# Missing value handling
ContrastsMissing$new(lfq_data, contrast_def)
```

**Key prophosqua functions:**
```r
# Data loading
data <- prophosqua::load_and_preprocess_data(dea_path)
data <- prophosqua::filter_contaminants(data)

# DPU calculation (protein-normalized PTM)
dpu_result <- prophosqua::test_diff(site_data, protein_data)

# Visualization
prophosqua::n_to_c_expression(data, contrast)
prophosqua::n_to_c_expression_multicontrast(data, contrasts)
```

### Data Manipulation

```r
library(dplyr)          # Data transformation
library(tidyr)          # Reshaping, pivoting
library(arrow)          # Parquet file I/O (fast, columnar storage)
library(readxl)         # Excel reading
library(writexl)        # Excel writing
library(yaml)           # Configuration files
```

### Enrichment Analysis

```r
# Ortholog mapping (for mouse studies)
library(biomaRt)        # Ensembl ortholog queries (cache results!)

# Gene set enrichment
library(fgsea)          # Fast GSEA implementation (PTM-SEA)
library(clusterProfiler) # GO, KEGG, MSigDB analysis
library(msigdbr)        # MSigDB gene sets

# Annotation databases
library(org.Mm.eg.db)   # Mouse
library(org.Hs.eg.db)   # Human

# Local KSEA
library(KSEAapp)        # KSEA.Scores() with NetworKIN database
```

### Visualization

```r
library(ggplot2)        # Base plotting
library(ggseqlogo)      # Sequence logos for kinase motifs
library(patchwork)      # Multi-panel composition
library(ComplexHeatmap) # Advanced heatmaps
library(circlize)       # Color palettes
library(enrichplot)     # Enrichment visualization
library(DT)             # Interactive tables
```

### Reporting

```r
library(rmarkdown)      # Document generation
library(bookdown)       # Extended Rmd (html_document2, pdf_document2)
library(knitr)          # Chunk options
library(here)           # Project-relative paths
```

## Analysis Workflow Architecture

### 1. Differential Expression Analysis (Foundation)

```
Input: Search engine output (Spectronaut BGS, FragPipe, MaxQuant)
       ↓
Normalization (VSN or median centering)
       ↓
Linear model fitting per site
       ↓
Contrast computation
       ↓
Output: Parquet (normalized abundances) + Excel (statistics)
```

### 2. DPA (Differential PTM Abundance)

Direct comparison of phosphosite signals between conditions.

**Key columns:**
- `diff.site`: log2 fold change
- `FDR.site`: FDR-adjusted p-value
- `gene_name.site`: Gene symbol

```r
# Join phospho to protein data
combined <- left_join(
  phospho_stats,
  protein_stats,
  by = c("protein_Id", "contrast"),
  suffix = c(".site", ".protein")
)
```

### 3. DPU (Differential PTM Usage)

Protein-normalized phosphorylation changes (stoichiometry).

**Calculation:**
```
diff_diff = log2(FC_PTM) - log2(FC_Protein)
SE_diff = sqrt(SE_PTM² + SE_Protein²)
```

**Key columns:**
- `diff_diff`: Protein-normalized log2 fold change
- `FDR_I`: FDR for interaction effect
- `protein_Id_site`: Site identifier

```r
dpu_result <- prophosqua::test_diff(
  site_data,
  protein_data,
  join_by = c("protein_Id", "contrast")
)
```

### 4. CorrectFirst (Alternative Approach)

Subtract protein abundance BEFORE statistical modeling:

```r
# On log scale
ptm_usage <- normalized_abundance.site - normalized_abundance.total

# Then fit linear model to corrected values
model <- build_model(ptm_usage_data, strategy)
```

### 5. KSEA Conversion

Convert results to KSEA format (Gene + Position):

```r
convert_to_ksea <- function(xlsx_file, sheet, output_dir,
                           diff_col, fdr_col, gene_col,
                           fdr_threshold = 0.05,
                           diff_threshold = 0.6) {
  data <- readxl::read_excel(xlsx_file, sheet = sheet)

  # Site-centric aggregation
  ksea_data <- data %>%
    filter(!is.na(!!sym(diff_col)), abs(!!sym(diff_col)) > diff_threshold) %>%
    group_by(Gene, Position) %>%
    summarize(
      FC = mean(!!sym(diff_col)),
      p_value = min(!!sym(fdr_col))
    )

  # Write per contrast
  for (contrast in unique(data$contrast)) {
    write_csv(
      filter(ksea_data, contrast == !!contrast),
      file.path(output_dir, paste0("KSEA_", contrast, ".csv"))
    )
  }
}
```

**Mouse to Human Mapping (cache results!):**

```r
get_mouse_human_orthologs <- function(cache_file, force_refresh = FALSE) {
  if (file.exists(cache_file) && !force_refresh) {
    return(readRDS(cache_file))
  }

  mart_mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  orthologs <- biomaRt::getBM(
    attributes = c("mgi_symbol", "hsapiens_homolog_associated_gene_name"),
    mart = mart_mouse
  )

  saveRDS(orthologs, cache_file)
  return(orthologs)
}
```

### 6. PTM-SEA Conversion (Recommended for Mouse)

Species-invariant using flanking sequences (±7 amino acids):

```r
convert_to_ptm_sea <- function(xlsx_file, sheet, output_dir,
                               diff_col, flanking_col = "flanking_sequence") {
  data <- readxl::read_excel(xlsx_file, sheet = sheet)

  # Create GCT format
  gct_data <- data %>%
    mutate(id = paste0(!!sym(flanking_col), "-p")) %>%
    select(id, Description = gene_name, value = !!sym(diff_col))

  write_gct(gct_data, file.path(output_dir, "ptmsea_input.gct"))
}

write_gct <- function(data, output_file, contrast_name) {
  # GCT v1.2 format
  con <- file(output_file, "w")
  writeLines("#1.2", con)
  writeLines(paste(nrow(data), 1, sep = "\t"), con)
  writeLines(paste("id", "Description", contrast_name, sep = "\t"), con)
  write.table(data, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(con)
}
```

## Snakemake Pipeline Pattern

```python
# Snakefile.smk
from datetime import date

DIR_OUT = f"PTM_{date.today()}"

CONTRASTS = [
    "treatment_vs_control",
    "condition_A_vs_B",
    "interaction_effect"
]

rule all:
    input:
        expand(f"{DIR_OUT}/PTM_DPA/KSEA/KSEA_{{contrast}}.csv", contrast=CONTRASTS),
        expand(f"{DIR_OUT}/PTM_DPA/PTMSEA/PTMSEA_{{contrast}}.gct", contrast=CONTRASTS),
        f"{DIR_OUT}/index.html"

rule dea:
    input:
        "Analysis_DEA.Rmd"
    output:
        f"{DIR_OUT}/PTM_DPA/Result_DPA.xlsx",
        f"{DIR_OUT}/PTM_DPU/Result_DPU.xlsx"
    shell:
        """
        Rscript -e "rmarkdown::render('Analysis_DEA.Rmd', params=list(output_dir='{DIR_OUT}'))"
        """

rule ksea_convert:
    input:
        f"{DIR_OUT}/PTM_DPA/Result_DPA.xlsx"
    output:
        f"{DIR_OUT}/PTM_DPA/KSEA/KSEA_{{contrast}}.csv"
    shell:
        """
        Rscript convert_to_KSEA.R --input {input} --contrast {wildcards.contrast} --output {output}
        """

rule gsea:
    input:
        f"{DIR_OUT}/PTM_DPA/DPA_GSEA_{{contrast}}.rnk"
    output:
        f"{DIR_OUT}/PTM_DPA/GSEA/GSEA_{{contrast}}.html"
    params:
        contrast = "{contrast}"
    shell:
        """
        Rscript -e "rmarkdown::render('GSEA_rnk.Rmd', params=list(rnk_file='{input}', contrast='{params.contrast}'), output_file='{output}')"
        """
```

## Column Name Mapping Reference

| Analysis | Diff Column | FDR Column | Gene Column | Site ID |
|----------|-------------|------------|-------------|---------|
| DPA | `diff.site` | `FDR.site` | `gene_name.site` | `protein_Id_site` |
| DPU | `diff_diff` | `FDR_I` | `gene_name.site` | `protein_Id_site` |
| CorrectFirst | `diff_diff` | `FDR_I` | `gene_name` | `protein_Id_site` |

**Handle via parameterization:**
```r
process_results <- function(data, diff_col, fdr_col, gene_col) {
  data %>%
    filter(!!sym(fdr_col) < 0.05) %>%
    arrange(desc(abs(!!sym(diff_col))))
}
```

## Output Directory Structure

```
PTM_YYYYMMDD/
├── PTM_DPA/
│   ├── Result_DPA.xlsx
│   ├── Analysis_DPA_n_to_c.html
│   ├── Analysis_DPA_seqlogo.html
│   ├── KSEA/
│   │   ├── KSEA_<contrast>.csv
│   │   └── *_ORA_SITE_*.txt
│   ├── PTMSEA/
│   │   └── PTMSEA_<contrast>.gct
│   └── GSEA/
│       ├── DPA_GSEA_<contrast>.rnk
│       └── GSEA_<contrast>.html
├── PTM_DPU/
│   └── (similar structure)
├── PTM_CF_DPU/
│   └── (CorrectFirst results)
└── index.html
```

## Best Practices

### 1. Cache Expensive Computations

```r
# Ortholog mapping
orthologs <- if (file.exists("cache/orthologs.rds")) {
  readRDS("cache/orthologs.rds")
} else {
  result <- fetch_orthologs()
  saveRDS(result, "cache/orthologs.rds")
  result
}
```

### 2. Validate Data Before Processing

```r
validate_phospho_data <- function(data) {
  required_cols <- c("protein_Id", "gene_name", "position", "flanking_sequence")
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  # Check flanking sequence format
  if (!all(nchar(data$flanking_sequence) == 15)) {
    warning("Some flanking sequences are not 15 characters")
  }
}
```

### 3. Site-Centric Aggregation

```r
# Multiple phosphopeptides may map to same site
aggregate_to_site <- function(data, diff_col, fdr_col) {
  data %>%
    group_by(protein_Id, gene_name, position) %>%
    summarize(
      diff = mean(!!sym(diff_col), na.rm = TRUE),
      fdr = min(!!sym(fdr_col), na.rm = TRUE),
      n_peptides = n()
    ) %>%
    ungroup()
}
```

### 4. Parameterized Rmd Documents

```yaml
# In Rmd YAML header
params:
  output_dir: "PTM_output"
  contrast: "treatment_vs_control"
  fdr_threshold: 0.05
  diff_threshold: 0.6
```

```r
# In Rmd body
results <- process_contrast(
  data,
  contrast = params$contrast,
  fdr_threshold = params$fdr_threshold
)
```

### 5. Generate ORA and GSEA Files

```r
generate_ora_gsea_files <- function(data, diff_col, fdr_col, gene_col,
                                    output_dir, contrast,
                                    fdr_threshold = 0.05,
                                    diff_threshold = 0.6) {
  sig_data <- data %>%
    filter(!!sym(fdr_col) < fdr_threshold)

  # ORA: up-regulated genes
  up_genes <- sig_data %>%
    filter(!!sym(diff_col) > diff_threshold) %>%
    pull(!!sym(gene_col)) %>%
    unique()
  writeLines(up_genes, file.path(output_dir, paste0(contrast, "_UP_ORA.txt")))

  # ORA: down-regulated genes
  down_genes <- sig_data %>%
    filter(!!sym(diff_col) < -diff_threshold) %>%
    pull(!!sym(gene_col)) %>%
    unique()
  writeLines(down_genes, file.path(output_dir, paste0(contrast, "_DOWN_ORA.txt")))

  # GSEA: ranked list
  rnk_data <- data %>%
    select(gene = !!sym(gene_col), score = !!sym(diff_col)) %>%
    group_by(gene) %>%
    summarize(score = mean(score, na.rm = TRUE)) %>%
    arrange(desc(score))
  write_tsv(rnk_data, file.path(output_dir, paste0(contrast, ".rnk")), col_names = FALSE)
}
```

## Sequence Logo Analysis

```r
library(ggseqlogo)

plot_kinase_motif <- function(sequences, title = "Phosphorylation Motif") {
  # sequences: character vector of ±7 flanking sequences (15 chars, center = phosphosite)

  # Build position frequency matrix
  pfm <- consensusMatrix(sequences)

  ggseqlogo(pfm, method = "prob") +
    labs(title = title) +
    theme_minimal()
}

plot_diff_logo <- function(up_sequences, down_sequences) {
  # Difference logo showing enriched/depleted amino acids
  pfm_up <- consensusMatrix(up_sequences)
  pfm_down <- consensusMatrix(down_sequences)

  diff_matrix <- pfm_up - pfm_down

  ggseqlogo(diff_matrix, method = "custom") +
    labs(title = "Upregulated - Downregulated Motif") +
    theme_minimal()
}
```

## Common Issues and Solutions

### 1. Missing Protein Data for DPU

```r
# Check match rate before DPU calculation
match_rate <- mean(phospho_data$protein_Id %in% protein_data$protein_Id)
message(sprintf("Phospho-protein match rate: %.1f%%", match_rate * 100))

if (match_rate < 0.5) {
  warning("Low match rate. Check protein ID format consistency.")
}
```

### 2. Mouse-Human Position Mismatch

KSEA databases use human positions. Mouse ortholog positions may differ.

**Solution:** Use PTM-SEA with flanking sequences (species-invariant).

### 3. Large Memory Usage

```r
# Use arrow/parquet for large datasets
library(arrow)

# Write intermediate results
write_parquet(large_data, "intermediate.parquet")

# Read only needed columns
data <- read_parquet("intermediate.parquet", col_select = c("protein_Id", "diff", "fdr"))
```

### 4. Contaminant Filtering

```r
filter_contaminants <- function(data) {
  data %>%
    filter(!grepl("^rev_|^CON_|^contaminant", protein_Id, ignore.case = TRUE))
}
```

## External Resources

- **KSEA App**: https://casecpb.shinyapps.io/ksea/
- **PTMsigDB**: https://github.com/broadinstitute/ssGSEA2.0
- **PhosphoSitePlus**: https://www.phosphosite.org/
- **MSigDB**: https://www.gsea-msigdb.org/
- **prolfqua documentation**: https://fgcz.github.io/prolfqua/
- **prophosqua documentation**: https://fgcz.github.io/prophosqua/

## Example Project Structure

```
project/
├── CLAUDE.md                    # Project documentation
├── Snakefile.smk               # Pipeline orchestration
├── config.yaml                 # Analysis parameters
├── Analysis_DEA.Rmd            # Main DEA workflow
├── Analysis_DPA_n_to_c.Rmd     # DPA visualization
├── Analysis_DPU_n_to_c.Rmd     # DPU visualization
├── Analysis_seqlogo.Rmd        # Sequence logo analysis
├── Analysis_KSEA.Rmd           # KSEA conversion
├── GSEA_rnk.Rmd               # GSEA analysis template
├── convert_to_KSEA.R          # Format conversion functions
├── SeqLogoDiff.R              # Visualization functions
├── DEA_input/                 # Raw data and config
│   ├── phospho_data.tsv
│   ├── protein_data.tsv
│   └── annotation.yaml
├── cache/                     # Cached computations
│   ├── mouse_human_orthologs.rds
│   └── ptmsigdb_v2.gmt
└── PTM_YYYYMMDD/              # Date-stamped outputs
```
