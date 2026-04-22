---
name: rmd-formatting
description: R Markdown formatting patterns: chunk display strategies (include=FALSE + inline R vs results='asis'), table rendering with knitr::kable, header hierarchy with tabsets, downloadable file embedding, and what to avoid (cat/print in visible chunks).
---

# R Markdown Report Formatting

## Core principle

Never use `cat()` or `print()` in visible chunks for user-facing output — it renders as monospace code blocks with `##` prefixes. Use `include=FALSE` + inline R, or `results='asis'`.

## Pattern 1: Simple statistics (preferred)

```r
# Hidden computation
```{r, include=FALSE}
n_cells <- ncol(obj)
n_genes <- nrow(obj)
```

**Dataset:** `r format(n_cells, big.mark = ",")` cells · `r format(n_genes, big.mark = ",")` genes
```

## Pattern 2: Tables

```r
```{r, include=FALSE}
sample_counts <- as.data.frame(table(obj$Sample))
colnames(sample_counts) <- c("Sample", "Cells")
```

```{r, echo=FALSE}
knitr::kable(sample_counts)
```
```

## Pattern 3: results='asis' (mixed compute + display in one chunk)

When a chunk must compute AND render text (e.g., QC loops):

**Required global setup:**
```r
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, comment = "")
```
```

Rules for `results='asis'` chunks:
1. `cat("text\n")` for plain text — renders as markdown
2. Tables: `knitr::kable(df, format = "markdown") |> cat(sep = "\n")` — NOT `print(df)`
3. End text blocks with `\n\n` for paragraph spacing

```r
```{r qc-summary, results='asis'}
n_pass <- sum(obj$pass_qc)
n_total <- ncol(obj)
cat(sprintf("**%d / %d cells pass QC (%.0f%%)**\n\n", n_pass, n_total, 100*n_pass/n_total))

tbl <- as.data.frame(table(obj$CellType, obj$pass_qc))
knitr::kable(tbl, format = "markdown") |> cat(sep = "\n")
```
```

## Pattern 4: Downloadable files

```r
# In a chunk with results='asis':
xfun::embed_file("path/to/file.xlsx", name = "output.xlsx", text = "Download Results (Excel)")
```

## Header hierarchy

- `# {.tabset}` — FIRST tabset only (level 1)
- `## Tab Name` — tabs within the first tabset; also all subsequent section headers
- `## {.tabset}` — any subsequent tabsets (level 2)
- `###` — subsections within tabs, or tabs of a level-2 tabset
- `####` — sub-subsections

**Start reports with a 1–2 sentence project description, then `# {.tabset}` directly. No "Executive Summary" or "Data Loading" sections at the top — data loading goes inside the first tab.**

```markdown
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, comment = "")
```

Brief project description here.

# {.tabset}

## Data Loading

## Analysis

### Subsection

## Session Information
```
```

## Quick reference

| Situation | Use |
|-----------|-----|
| Simple stats (1–3 values) | `include=FALSE` + inline R |
| Complex/conditional text | `results='asis'` + `cat()` |
| Tables (always) | `knitr::kable()` |
| Never | `print()`, `print(table())` in visible chunks |
