---
name: marimo-compact
description: Compact marimo notebook patterns. Reactive cells, UI elements, data visualization.
---

# Marimo Essentials

## Cell Structure
```python
@app.cell
def _():
    import marimo as mo
    return

@app.cell
def _():
    slider = mo.ui.slider(0, 100, value=50, label="Value")
    slider  # Display it
    return

@app.cell
def _():
    result = slider.value * 2  # Reactive: re-runs when slider changes
    mo.md(f"Result: {result}")
    return
```

## Reactivity Rules
1. Variables auto-propagate: changing `x` re-runs cells using `x`
2. No redeclaration: each variable defined once across all cells
3. UI `.value` in separate cell from definition
4. `_var` = cell-local variable
5. Last expression displays automatically

## Common UI Elements
```python
mo.ui.slider(start, stop, value=None, label=None)
mo.ui.dropdown(options, value=None, label=None)
mo.ui.text(value="", label=None)
mo.ui.checkbox(label="", value=False)
mo.ui.button(value=None, kind="primary")
mo.ui.table(data, sortable=True, filterable=True)
mo.ui.dataframe(df)  # Editable
mo.ui.data_explorer(df)  # Interactive explorer
mo.ui.form(element, label="")  # Wrap in form
```

## Layouts
```python
mo.md("# Title")  # Markdown
mo.hstack([a, b, c])  # Horizontal
mo.vstack([top, mid, bot])  # Vertical
mo.tabs({"Tab 1": content1, "Tab 2": content2})
mo.accordion({"Section": content})
```

## Data & SQL
```python
# DataFrames display automatically
import polars as pl
df = pl.read_csv("data.csv")
df  # Shows as table

# SQL with DuckDB
result = mo.sql(f"SELECT * FROM df WHERE x > 100")
```

## Visualization
- **matplotlib**: `plt.gca()` as last expression (not `plt.show()`)
- **plotly**: return figure object directly
- **altair**: return chart object; accepts polars DataFrames

## Control Flow
```python
mo.stop(condition, mo.md("Stopped"))  # Stop cell execution
mo.output.replace(content)  # Replace output
```

## Pitfalls to Avoid
- Circular dependencies between cells
- Accessing `.value` in same cell as UI definition
- Using `global` keyword
- Not returning visualization object as last expression

## CLI
```bash
marimo edit notebook.py  # Edit
marimo run notebook.py   # Run as app
marimo check --fix       # Fix formatting issues
```
