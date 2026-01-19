---
name: bookdown
description: Cross-referencing in R Markdown documents using bookdown. Use when creating Rmd files that need references to figures, tables, equations, or sections. Triggers on requests for academic reports, technical documents, or any R Markdown with numbered cross-references.
---

# Bookdown Cross-Referencing

## Output Formats

Use `bookdown::` output formats for cross-referencing support:

```yaml
output:
  bookdown::html_document2:
    toc: true
  bookdown::pdf_document2:
    toc: true
```

Other formats: `gitbook`, `word_document2`, `epub_book`.

## Figures

Chunk name becomes the reference ID. Requires `fig.cap`:

```r
```{r myplot, fig.cap="Distribution of values."}
plot(x)
```
```

Reference: `Figure \@ref(fig:myplot)` → "Figure 1"

### Long Captions (Text Reference)

Define caption outside chunk, reference inside:

```markdown
(ref:myplot) This is a very long caption with **markdown** and citations [@smith2020].

```{r myplot, fig.cap="(ref:myplot)"}
plot(x)
```
```

## Tables

Use `knitr::kable()` with `caption`:

```r
```{r mytable}
knitr::kable(df, caption = "Summary statistics.")
```
```

Reference: `Table \@ref(tab:mytable)` → "Table 1"

## Sections

Headers auto-generate IDs from text (lowercase, hyphens):

```markdown
# My Analysis Section
See Section \@ref(my-analysis-section).
```

Custom ID: `# My Section {#custom-id}`

## Equations

```markdown
\begin{equation}
y = mx + b
(\#eq:linear)
\end{equation}

See Equation \@ref(eq:linear).
```

## Quick Reference

| Element  | Chunk/Label | Reference Syntax |
|----------|-------------|------------------|
| Figure   | `myplot`    | `\@ref(fig:myplot)` |
| Table    | `mytable`   | `\@ref(tab:mytable)` |
| Section  | `my-section`| `\@ref(my-section)` |
| Equation | `eq:model`  | `\@ref(eq:model)` |

## Common Issues

- **References show as `??`**: Missing `fig.cap` or wrong output format
- **Duplicate labels**: Each chunk name must be unique
- **Backslash required**: Use `\@ref()` not `@ref()`
