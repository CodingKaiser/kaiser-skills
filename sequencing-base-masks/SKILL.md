---
name: sequencing-base-masks
description: Sequencing read requirements and bcl2fastq base masks for 10x Genomics libraries: chemistry-specific BasesMaskCR values, ATAC I2 barcode requirement (capital I not lowercase y), Visium FFPE vs frozen R2 length difference.
---

# Sequencing Read Requirements & Base Masks

## Key concept: `I` vs `y`

- **`I` (capital)** = output as index FASTQ (I1/I2) — used for sample demultiplexing
- **`y` (lowercase)** = output as regular read FASTQ (R1, R2, R3…)

For 10x ATAC, the i5 read IS the cell barcode — **must use capital `I`** to produce I2 FASTQ, or CellRanger ARC fails.

## Base Masks by Chemistry

### 3' / 5' Gene Expression

| Chemistry | BasesMaskCR |
|-----------|-------------|
| 3' v2 (legacy) | `y26n*,I8n*,y*` |
| 3' v3/v4/GEM-X | `y28n*,I10n*,I10n*,y*` |
| 5' v2/HT (legacy) | `y26n*,I10n*,I10n*,y*` |
| 5' v3 / Universal 5' | `y28n*,I10n*,I10n*,y*` |
| 5' R1-only | `y41n*,I10n*,I10n*` |
| 5' R1-only v3 | `y43n*,I10n*,I10n*` |

### VDJ / Feature Barcoding / CellPlex

| Chemistry | BasesMaskCR |
|-----------|-------------|
| VDJ v2/v3 | `y28n*,I10n*,I10n*,y91n*` |
| Feature Barcode (CITE-seq) | `y28n*,I10n*,I10n*,y*` |
| CellPlex / CMO | `y28n*,I10n*,I10n*,y*` |

### ATAC & Multiome

| Chemistry | BasesMaskCR | Notes |
|-----------|-------------|-------|
| ATAC solo v1/v2 | `y50n*,I8n*,I16n*,y50n*` | i5=16bp barcode → capital I |
| Multiome GEX | `y28n*,I10n*,I10n*,y*` | |
| Multiome ATAC | `y50n*,I8n*,I24n*,y49n*` | i5=24bp (16bp BC + 8bp linker) → capital I |

**Do NOT pool ATAC libraries with other 10x types** — incompatible i5.

### Spatial (Visium)

| Chemistry | BasesMaskCR | Notes |
|-----------|-------------|-------|
| Visium frozen | `y28n*,I10n*,I10n*,y90n*` | Kit TT |
| Visium FFPE | `y28n*,I10n*,I10n*,y50n*` | R2=50bp — shorter insert |
| Visium CytAssist | `y28n*,I10n*,I10n*,y90n*` | Same as frozen (NOT FFPE) |
| Visium HD | `y43n*,I10n*,I10n*,y151n*` | Extended R1 |

### Fixed RNA (Flex)

| Chemistry | BasesMaskCR |
|-----------|-------------|
| Flex GE v2 standard | `y28n*,I10n*,I10n*,y90n*` |
| Flex GE v2 R1-long | `y54n*,I10n*,I10n*,y50n*` |
| Flex Protein Expression | `y28n*,I10n*,I10n*,y90n*` |

## Recommended Sequencing Depths

| Chemistry | Depth |
|-----------|-------|
| 3' / 5' GEX, Multiome GEX | ≥20,000 pairs/cell |
| ATAC solo / Multiome ATAC | ≥25,000 pairs/nucleus |
| Visium frozen / CytAssist | ≥50,000 pairs/spot |
| Visium FFPE | ≥25,000 pairs/spot |
| Flex GE | ≥5,000–10,000 pairs/cell |
