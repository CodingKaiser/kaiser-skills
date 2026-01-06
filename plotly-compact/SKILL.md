---
name: plotly-compact
description: Compact Plotly visualization patterns. Express for quick plots, Graph Objects for control.
---

# Plotly Essentials

## Express vs Graph Objects
- **Express**: Quick exploratory plots from DataFrames
- **Graph Objects**: Fine-grained control, complex multi-trace figures

## Common Charts (Express)
```python
import plotly.express as px

px.scatter(df, x="x", y="y", color="cat", size="val", hover_data=["name"])
px.line(df, x="date", y="value", color="series")
px.bar(df, x="cat", y="val", color="group", barmode="group")  # or "stack"
px.histogram(df, x="val", color="group", barmode="overlay", opacity=0.7)
px.box(df, x="group", y="val", points="all")
px.violin(df, x="group", y="val", box=True)
px.imshow(corr_matrix, text_auto=".2f", color_continuous_scale="RdBu_r")
px.scatter_matrix(df, dimensions=["a", "b", "c"], color="cat")
```

## Graph Objects (when needed)
```python
import plotly.graph_objects as go
from plotly.subplots import make_subplots

fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name="series"))
fig.update_layout(title="Title", template="plotly_white")
fig.update_traces(marker=dict(size=10))
```

## Subplots
```python
fig = make_subplots(rows=2, cols=2, subplot_titles=["A", "B", "C", "D"])
fig.add_trace(go.Scatter(x=x, y=y), row=1, col=1)

# Secondary y-axis
fig = make_subplots(specs=[[{"secondary_y": True}]])
fig.add_trace(go.Scatter(...), secondary_y=False)
fig.add_trace(go.Bar(...), secondary_y=True)
```

## Customization
```python
fig.update_layout(
    template="plotly_white",  # or plotly_dark, ggplot2, seaborn
    legend=dict(orientation="h", yanchor="bottom", y=1.02)
)
fig.add_hline(y=threshold, line_dash="dash")
fig.add_annotation(x=2, y=5, text="Note", showarrow=True)
fig.update_xaxes(type="log", title="X Label")
```

## Export
```python
fig.write_html("plot.html")
fig.write_image("plot.png", scale=2)  # requires kaleido
fig.write_image("plot.svg")
```

## Network Graphs
```python
import networkx as nx

G = nx.karate_club_graph()
pos = nx.spring_layout(G, seed=42)

# Edges
edge_x, edge_y = [], []
for u, v in G.edges():
    x0, y0 = pos[u]; x1, y1 = pos[v]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])

edge_trace = go.Scatter(x=edge_x, y=edge_y, mode="lines", line=dict(width=0.5, color="#888"))
node_trace = go.Scatter(x=[pos[n][0] for n in G.nodes()], y=[pos[n][1] for n in G.nodes()],
                        mode="markers", marker=dict(size=10, color=list(dict(G.degree()).values())))

fig = go.Figure([edge_trace, node_trace])
fig.update_layout(showlegend=False, xaxis=dict(showgrid=False, showticklabels=False),
                  yaxis=dict(showgrid=False, showticklabels=False))
```

## Tips
- Large data: `px.scatter(..., render_mode="webgl")`
- Colorblind-safe: `color_discrete_sequence=px.colors.qualitative.Safe`
