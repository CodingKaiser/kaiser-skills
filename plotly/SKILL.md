---
name: plotly
description: Plotly visualization patterns for statistical and scientific charts. Use when creating interactive visualizations, statistical plots (scatter, box, violin, heatmaps), UpSet plots for set intersections, network graphs, or exporting figures to HTML/PNG/PDF/SVG formats. Covers both Plotly Express (high-level) and Graph Objects (low-level) APIs.
---

# Plotly Visualization Guide

Comprehensive patterns for creating statistical and scientific visualizations with Plotly, covering both the high-level Plotly Express API and the low-level Graph Objects API.

## Quick Reference

```python
# Essential imports
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Quick Express charts
fig = px.scatter(df, x="x", y="y", color="category")
fig = px.histogram(df, x="values", nbins=30)
fig = px.box(df, x="group", y="value")

# Show and export
fig.show()                           # Interactive display
fig.write_html("plot.html")          # Interactive HTML
fig.write_image("plot.png")          # Static image (requires kaleido)
```

## When to Use Express vs Graph Objects

| Use Case | Recommended API |
|----------|----------------|
| Quick exploratory plots | Express |
| Standard chart types with DataFrames | Express |
| Complex multi-trace figures | Graph Objects |
| Fine-grained customization | Graph Objects |
| Animations and transitions | Express (simpler) or GO |
| Custom trace types | Graph Objects |

**Rule of thumb:** Start with Express. Switch to Graph Objects when you need more control.

## Core Concepts

### Figure Structure

Every Plotly figure has two main components:

```python
fig = go.Figure(
    data=[...],     # List of traces (the actual chart data)
    layout={...}    # Layout configuration (axes, title, legend, etc.)
)
```

### Traces

Traces are the individual data series in a chart:

```python
# Graph Objects trace
trace = go.Scatter(x=[1, 2, 3], y=[4, 5, 6], mode="lines+markers")

# Express automatically creates traces from DataFrame columns
fig = px.scatter(df, x="x", y="y", color="category")  # One trace per category
```

### Updating Figures

```python
# Update layout
fig.update_layout(
    title="My Chart",
    xaxis_title="X Axis",
    yaxis_title="Y Axis",
    template="plotly_white"
)

# Update all traces
fig.update_traces(marker=dict(size=10))

# Chain updates
fig.update_layout(...).update_traces(...).update_xaxes(...)
```

## Statistical Charts

### Scatter Plots

```python
# Basic scatter with Express
fig = px.scatter(df, x="x", y="y")

# With color, size, and hover data
fig = px.scatter(
    df,
    x="x",
    y="y",
    color="category",
    size="magnitude",
    hover_data=["name", "value"],
    title="Scatter Plot"
)

# Scatter matrix (pairplot)
fig = px.scatter_matrix(
    df,
    dimensions=["col1", "col2", "col3"],
    color="category"
)

# Graph Objects for full control
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=df["x"],
    y=df["y"],
    mode="markers",
    marker=dict(
        size=10,
        color=df["value"],
        colorscale="Viridis",
        showscale=True
    ),
    text=df["label"],
    hovertemplate="<b>%{text}</b><br>x: %{x}<br>y: %{y}<extra></extra>"
))
```

### Line Charts

```python
# Basic line chart
fig = px.line(df, x="date", y="value")

# Multiple lines with grouping
fig = px.line(df, x="date", y="value", color="series", line_dash="type")

# With confidence intervals using Graph Objects
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=df["x"], y=df["upper"],
    mode="lines", line=dict(width=0),
    showlegend=False
))
fig.add_trace(go.Scatter(
    x=df["x"], y=df["lower"],
    mode="lines", line=dict(width=0),
    fill="tonexty", fillcolor="rgba(68, 68, 68, 0.3)",
    showlegend=False
))
fig.add_trace(go.Scatter(
    x=df["x"], y=df["mean"],
    mode="lines", name="Mean"
))
```

### Bar Charts

```python
# Vertical bar chart
fig = px.bar(df, x="category", y="value")

# Grouped bars
fig = px.bar(df, x="category", y="value", color="group", barmode="group")

# Stacked bars
fig = px.bar(df, x="category", y="value", color="group", barmode="stack")

# Horizontal bar chart
fig = px.bar(df, x="value", y="category", orientation="h")

# With error bars
fig = px.bar(df, x="category", y="mean", error_y="std")
```

### Histograms

```python
# Basic histogram
fig = px.histogram(df, x="value", nbins=30)

# Overlaid histograms
fig = px.histogram(
    df, x="value", color="group",
    barmode="overlay", opacity=0.7
)

# Normalized histogram (density)
fig = px.histogram(df, x="value", histnorm="probability density")

# 2D histogram (heatmap)
fig = px.density_heatmap(df, x="x", y="y", nbinsx=30, nbinsy=30)
```

### Box Plots

```python
# Basic box plot
fig = px.box(df, x="group", y="value")

# With individual points
fig = px.box(df, x="group", y="value", points="all")

# Notched box plot (shows confidence interval for median)
fig = px.box(df, x="group", y="value", notched=True)

# Grouped box plots
fig = px.box(df, x="group", y="value", color="subgroup")
```

### Violin Plots

```python
# Basic violin plot
fig = px.violin(df, x="group", y="value")

# With box plot inside
fig = px.violin(df, x="group", y="value", box=True)

# With individual points
fig = px.violin(df, x="group", y="value", points="all")

# Split violin (comparing two groups)
fig = px.violin(
    df, x="category", y="value", color="group",
    violinmode="overlay"
)
```

### Heatmaps and Correlation Matrices

```python
# Basic heatmap
fig = px.imshow(matrix, text_auto=True)

# Correlation matrix
corr = df.corr()
fig = px.imshow(
    corr,
    text_auto=".2f",
    color_continuous_scale="RdBu_r",
    zmin=-1, zmax=1,
    title="Correlation Matrix"
)

# Clustered heatmap with Graph Objects
import scipy.cluster.hierarchy as sch

# Compute hierarchical clustering
linkage = sch.linkage(matrix, method="ward")
order = sch.leaves_list(linkage)

fig = go.Figure(go.Heatmap(
    z=matrix[order][:, order],
    x=labels[order],
    y=labels[order],
    colorscale="Viridis"
))
```

### Strip and Swarm Plots

```python
# Strip plot (jittered points)
fig = px.strip(df, x="group", y="value", color="group")

# Combined with box plot
fig = px.box(df, x="group", y="value")
fig.add_trace(go.Scatter(
    x=df["group"],
    y=df["value"],
    mode="markers",
    marker=dict(size=5, opacity=0.5),
    showlegend=False
))
```

## UpSet Plots for Set Intersections

UpSet plots visualize set intersections more effectively than Venn diagrams for more than 3 sets.

### Using upsetplot Library with Plotly

```python
from upsetplot import from_memberships, UpSet
import matplotlib.pyplot as plt

# Prepare data as membership lists
data = from_memberships([
    ["Set A"],
    ["Set A", "Set B"],
    ["Set B", "Set C"],
    ["Set A", "Set B", "Set C"],
], data=[10, 5, 3, 2])

# Create matplotlib UpSet, then convert if needed
upset = UpSet(data, show_counts=True)
upset.plot()
```

### Manual UpSet Plot with Plotly

```python
def create_upset_plot(sets_dict: dict[str, set]) -> go.Figure:
    """
    Create an UpSet plot from a dictionary of sets.

    Args:
        sets_dict: Dictionary mapping set names to sets of elements
    """
    from itertools import combinations

    set_names = list(sets_dict.keys())
    n_sets = len(set_names)

    # Calculate all intersections
    intersections = []
    for r in range(1, n_sets + 1):
        for combo in combinations(range(n_sets), r):
            # Elements in all sets of this combination
            intersection = set.intersection(*[sets_dict[set_names[i]] for i in combo])
            # Exclude elements in other sets (exclusive intersection)
            for i in range(n_sets):
                if i not in combo:
                    intersection = intersection - sets_dict[set_names[i]]
            if intersection:
                intersections.append({
                    "sets": combo,
                    "count": len(intersection),
                    "label": " & ".join([set_names[i] for i in combo])
                })

    # Sort by count descending
    intersections.sort(key=lambda x: x["count"], reverse=True)

    # Create subplot figure
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.7, 0.3],
        vertical_spacing=0.02,
        shared_xaxes=True
    )

    # Top: bar chart of intersection sizes
    fig.add_trace(go.Bar(
        x=list(range(len(intersections))),
        y=[i["count"] for i in intersections],
        marker_color="steelblue",
        showlegend=False
    ), row=1, col=1)

    # Bottom: dot matrix showing set membership
    for idx, inter in enumerate(intersections):
        for set_idx in range(n_sets):
            fig.add_trace(go.Scatter(
                x=[idx],
                y=[set_idx],
                mode="markers",
                marker=dict(
                    size=12,
                    color="black" if set_idx in inter["sets"] else "lightgray"
                ),
                showlegend=False
            ), row=2, col=1)

        # Connect dots for multi-set intersections
        if len(inter["sets"]) > 1:
            fig.add_trace(go.Scatter(
                x=[idx, idx],
                y=[min(inter["sets"]), max(inter["sets"])],
                mode="lines",
                line=dict(color="black", width=2),
                showlegend=False
            ), row=2, col=1)

    # Update layout
    fig.update_layout(
        title="UpSet Plot",
        height=500,
        yaxis2=dict(
            ticktext=set_names,
            tickvals=list(range(n_sets)),
            tickmode="array"
        ),
        xaxis2=dict(showticklabels=False)
    )

    return fig

# Usage
sets = {
    "Set A": {1, 2, 3, 4, 5},
    "Set B": {3, 4, 5, 6, 7},
    "Set C": {5, 6, 7, 8, 9}
}
fig = create_upset_plot(sets)
fig.show()
```

## Network Graphs

Network graphs visualize relationships between entities (nodes connected by edges). Plotly renders networks using scatter traces for nodes and line traces for edges, combined with NetworkX for graph data structures and layout algorithms.

### When to Use Which Tool

| Need | Recommended Tool |
|------|------------------|
| Integration with other Plotly charts | Plotly + NetworkX |
| Static publication figures | Plotly + NetworkX |
| Quick exploration in Jupyter | [gravis](https://robert-haas.github.io/gravis-docs/) |
| Dash web applications | [Dash Cytoscape](https://dash.plotly.com/cytoscape) |

### Basic Network Graph

```python
import plotly.graph_objects as go
import networkx as nx

# Create graph
G = nx.karate_club_graph()

# Compute layout
pos = nx.spring_layout(G, seed=42)

# Build edge traces (use None to create disconnected line segments)
edge_x, edge_y = [], []
for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    mode="lines",
    line=dict(width=0.5, color="#888"),
    hoverinfo="none"
)

# Build node trace
node_x = [pos[node][0] for node in G.nodes()]
node_y = [pos[node][1] for node in G.nodes()]

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode="markers",
    hoverinfo="text",
    marker=dict(
        size=10,
        color="steelblue",
        line=dict(width=1, color="white")
    )
)

# Create figure with clean layout (no axes)
fig = go.Figure(
    data=[edge_trace, node_trace],
    layout=go.Layout(
        title="Network Graph",
        showlegend=False,
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        template="plotly_white"
    )
)
fig.show()
```

### Layout Algorithms

NetworkX provides several layout algorithms:

```python
# Force-directed (spring) layout - good general purpose
pos = nx.spring_layout(G, seed=42, k=1/sqrt(len(G.nodes())))

# Kamada-Kawai - minimizes edge length differences
pos = nx.kamada_kawai_layout(G)

# Circular - nodes arranged in a circle
pos = nx.circular_layout(G)

# Shell - concentric circles by node groups
pos = nx.shell_layout(G, nlist=[center_nodes, outer_nodes])

# Spectral - uses graph eigenvalues
pos = nx.spectral_layout(G)
```

### Node Styling by Metrics

```python
# Color nodes by degree (number of connections)
degrees = dict(G.degree())
node_colors = [degrees[node] for node in G.nodes()]

# Size nodes by centrality
centrality = nx.betweenness_centrality(G)
node_sizes = [20 + 50 * centrality[node] for node in G.nodes()]

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode="markers",
    marker=dict(
        size=node_sizes,
        color=node_colors,
        colorscale="Viridis",
        showscale=True,
        colorbar=dict(title="Degree"),
        line=dict(width=1, color="white")
    ),
    text=[f"Node {n}<br>Degree: {degrees[n]}" for n in G.nodes()],
    hoverinfo="text"
)
```

### Bipartite Networks

For networks with two distinct node types (e.g., terms and proteins):

```python
def bipartite_layout(G: nx.Graph, left_nodes: list, right_nodes: list) -> dict:
    """Create a bipartite layout with left and right columns."""
    pos = {}
    for i, node in enumerate(left_nodes):
        pos[node] = (-1, i / len(left_nodes))
    for i, node in enumerate(right_nodes):
        pos[node] = (1, i / len(right_nodes))
    return pos

# Style different node types differently
left_trace = go.Scatter(
    x=[pos[n][0] for n in left_nodes],
    y=[pos[n][1] for n in left_nodes],
    mode="markers+text",
    marker=dict(size=15, color="steelblue"),
    text=left_nodes,
    textposition="middle left"
)

right_trace = go.Scatter(
    x=[pos[n][0] for n in right_nodes],
    y=[pos[n][1] for n in right_nodes],
    mode="markers",
    marker=dict(size=8, color="coral")
)
```

### Node Border Colors for Categories

Use marker line color to encode a second categorical variable:

```python
# Map categories to border colors
border_colors = {
    "up": "red",
    "down": "blue",
    "neutral": "gray"
}

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode="markers",
    marker=dict(
        size=12,
        color=node_fill_colors,      # Fill by continuous value
        colorscale="RdYlGn",
        line=dict(
            width=2,
            color=[border_colors[cat] for cat in node_categories]
        )
    )
)
```

### Interactive Hover Templates

```python
# Rich hover information
node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode="markers",
    marker=dict(size=10, color="steelblue"),
    customdata=[[degrees[n], centrality[n]] for n in G.nodes()],
    hovertemplate=(
        "<b>%{text}</b><br>"
        "Degree: %{customdata[0]}<br>"
        "Centrality: %{customdata[1]:.3f}"
        "<extra></extra>"  # Hides secondary box
    ),
    text=list(G.nodes())
)
```

### Large Network Performance

For networks with >1000 nodes:

```python
# Reduce visual complexity
edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    mode="lines",
    line=dict(width=0.3, color="rgba(150,150,150,0.3)"),
    hoverinfo="none"
)

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode="markers",
    marker=dict(size=3),  # Smaller markers
    hoverinfo="none"      # Disable hover for performance
)

# Consider filtering to show only important nodes/edges
important_nodes = [n for n in G.nodes() if degrees[n] > threshold]
subgraph = G.subgraph(important_nodes)
```

### Complete Network Visualization Function

```python
def plot_network(
    G: nx.Graph,
    layout: str = "spring",
    node_color: str | list = "steelblue",
    node_size: int | list = 10,
    title: str = "Network Graph"
) -> go.Figure:
    """Create an interactive network visualization."""
    # Compute layout
    layout_funcs = {
        "spring": nx.spring_layout,
        "kamada_kawai": nx.kamada_kawai_layout,
        "circular": nx.circular_layout,
    }
    pos = layout_funcs.get(layout, nx.spring_layout)(G)

    # Build edges
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line=dict(width=0.5, color="#888"),
        hoverinfo="none"
    )

    # Build nodes
    node_x = [pos[n][0] for n in G.nodes()]
    node_y = [pos[n][1] for n in G.nodes()]

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode="markers",
        marker=dict(
            size=node_size,
            color=node_color,
            line=dict(width=1, color="white")
        ),
        text=list(G.nodes()),
        hoverinfo="text"
    )

    return go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            title=title,
            showlegend=False,
            hovermode="closest",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            template="plotly_white"
        )
    )

# Usage
fig = plot_network(G, layout="kamada_kawai", node_color=list(degrees.values()))
fig.show()
```

## Subplots and Multiple Axes

### Creating Subplots

```python
from plotly.subplots import make_subplots

# Basic grid
fig = make_subplots(rows=2, cols=2, subplot_titles=["A", "B", "C", "D"])

fig.add_trace(go.Scatter(x=[1, 2], y=[1, 2]), row=1, col=1)
fig.add_trace(go.Bar(x=[1, 2], y=[3, 4]), row=1, col=2)
fig.add_trace(go.Scatter(x=[1, 2], y=[2, 1]), row=2, col=1)
fig.add_trace(go.Histogram(x=[1, 1, 2, 3, 3, 3]), row=2, col=2)

# Different sized subplots
fig = make_subplots(
    rows=2, cols=2,
    column_widths=[0.7, 0.3],
    row_heights=[0.3, 0.7],
    specs=[[{"colspan": 2}, None],
           [{}, {}]]
)

# Shared axes
fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
```

### Secondary Y-Axis

```python
fig = make_subplots(specs=[[{"secondary_y": True}]])

fig.add_trace(
    go.Scatter(x=[1, 2, 3], y=[40, 50, 60], name="Series A"),
    secondary_y=False
)
fig.add_trace(
    go.Bar(x=[1, 2, 3], y=[4, 5, 6], name="Series B"),
    secondary_y=True
)

fig.update_yaxes(title_text="Series A", secondary_y=False)
fig.update_yaxes(title_text="Series B", secondary_y=True)
```

## Customization Patterns

### Themes and Templates

```python
# Built-in templates
fig.update_layout(template="plotly_white")   # Clean white background
fig.update_layout(template="plotly_dark")    # Dark theme
fig.update_layout(template="ggplot2")        # ggplot2 style
fig.update_layout(template="seaborn")        # Seaborn style

# Set default template globally
import plotly.io as pio
pio.templates.default = "plotly_white"
```

### Colors and Color Scales

```python
# Named color sequences
fig = px.scatter(df, x="x", y="y", color="category",
                 color_discrete_sequence=px.colors.qualitative.Set2)

# Continuous color scales
fig = px.scatter(df, x="x", y="y", color="value",
                 color_continuous_scale="Viridis")

# Common scales: "Viridis", "Plasma", "Inferno", "RdBu", "Blues", "Reds"

# Custom discrete colors
color_map = {"A": "#1f77b4", "B": "#ff7f0e", "C": "#2ca02c"}
fig = px.scatter(df, x="x", y="y", color="category", color_discrete_map=color_map)
```

### Annotations and Shapes

```python
# Text annotations
fig.add_annotation(
    x=2, y=5,
    text="Important point",
    showarrow=True,
    arrowhead=2
)

# Horizontal/vertical lines
fig.add_hline(y=3, line_dash="dash", line_color="red")
fig.add_vline(x=5, line_dash="dot", annotation_text="Threshold")

# Rectangular regions
fig.add_vrect(x0=2, x1=4, fillcolor="green", opacity=0.2, line_width=0)
fig.add_hrect(y0=1, y1=3, fillcolor="blue", opacity=0.1, line_width=0)

# Shapes
fig.add_shape(
    type="circle",
    x0=1, y0=1, x1=3, y1=3,
    line=dict(color="red", dash="dash")
)
```

### Legend Configuration

```python
fig.update_layout(
    legend=dict(
        title="Categories",
        orientation="h",        # Horizontal legend
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1
    )
)

# Hide legend
fig.update_layout(showlegend=False)

# Legend inside plot
fig.update_layout(legend=dict(x=0.02, y=0.98))
```

### Axis Formatting

```python
# Log scale
fig.update_xaxes(type="log")
fig.update_yaxes(type="log")

# Date axis formatting
fig.update_xaxes(
    dtick="M1",              # Monthly ticks
    tickformat="%b %Y"       # Format: "Jan 2024"
)

# Reversed axis
fig.update_yaxes(autorange="reversed")

# Axis range
fig.update_xaxes(range=[0, 100])

# Grid styling
fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="lightgray")
```

## Export Workflows

### Interactive HTML

```python
# Full HTML file (standalone)
fig.write_html("plot.html", include_plotlyjs=True)

# Smaller HTML (requires CDN)
fig.write_html("plot.html", include_plotlyjs="cdn")

# HTML div only (for embedding)
html_div = fig.to_html(full_html=False, include_plotlyjs=False)
```

### Static Images

Requires `kaleido` package: `pip install kaleido`

```python
# PNG (raster)
fig.write_image("plot.png", scale=2)  # scale=2 for higher resolution

# SVG (vector)
fig.write_image("plot.svg")

# PDF (vector)
fig.write_image("plot.pdf")

# Specify dimensions
fig.write_image("plot.png", width=800, height=600, scale=2)
```

### Batch Export

```python
def export_figure(fig: go.Figure, basename: str, formats: list[str] = None):
    """Export figure to multiple formats."""
    if formats is None:
        formats = ["html", "png", "svg"]

    for fmt in formats:
        if fmt == "html":
            fig.write_html(f"{basename}.html")
        else:
            fig.write_image(f"{basename}.{fmt}", scale=2)
```

## Best Practices

### General Guidelines

- [ ] Start with Express for quick exploration, switch to Graph Objects for customization
- [ ] Use `template="plotly_white"` for publication-ready figures
- [ ] Always add axis labels and titles
- [ ] Use `hover_data` to add context without cluttering the plot
- [ ] Use `color_discrete_sequence` for consistent categorical colors
- [ ] Export to SVG/PDF for publications, PNG for presentations

### Performance Tips

- [ ] For large datasets (>10k points), use `px.scatter` with `render_mode="webgl"`
- [ ] Use `fig.update_traces(hoverinfo="skip")` to disable hover on dense plots
- [ ] Consider downsampling or aggregating before plotting millions of points

### Accessibility

- [ ] Use colorblind-friendly palettes: `px.colors.qualitative.Safe`
- [ ] Ensure sufficient contrast between colors
- [ ] Add text labels or patterns in addition to color encoding
- [ ] Include descriptive titles and axis labels

## Resources

- [Plotly Python Documentation](https://plotly.com/python/)
- [Plotly Express API Reference](https://plotly.com/python-api-reference/plotly.express.html)
- [Graph Objects Reference](https://plotly.com/python-api-reference/plotly.graph_objects.html)
- [Plotly Color Scales](https://plotly.com/python/builtin-colorscales/)
- [Plotly Network Graphs](https://plotly.com/python/network-graphs/)
- [NetworkX Documentation](https://networkx.org/documentation/stable/)
- [gravis](https://robert-haas.github.io/gravis-docs/) (interactive graph visualization)
- [Dash Cytoscape](https://dash.plotly.com/cytoscape) (for Dash applications)
- [UpSetPlot Library](https://upsetplot.readthedocs.io/) (for advanced UpSet plots)

## Common Issues

### "kaleido" not found for image export

```bash
pip install -U kaleido
```

### Figure not displaying in Jupyter

```python
# Ensure notebook renderer is set
import plotly.io as pio
pio.renderers.default = "notebook"  # or "jupyterlab", "vscode"
```

### Slow rendering with large datasets

```python
# Use WebGL renderer for scatter plots
fig = px.scatter(df, x="x", y="y", render_mode="webgl")

# Or reduce data points
fig = px.scatter(df.sample(10000), x="x", y="y")
```
