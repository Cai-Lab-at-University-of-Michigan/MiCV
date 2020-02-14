import dash_core_components as dcc
import plotly.graph_objs as go

import pandas as pd

from plotting.plotting_parameters import scale

### PROCESSING TAB CUSTOM COMPONENTS ###
def plot_pseudotime_UMAP():
    g = dcc.Graph(
        id='pseudotime_calculation_UMAP_plot',
        figure={
            'data': [
                go.Scattergl(
                    x=None,
                    y=None,
                    text=str("NULL"),
                    mode='markers',
                    opacity=0.7,
                    marker={
                        "line_width": 1
                    },
                    name=str("NULL")
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'UMAP1'},
                yaxis={'title': "UMAP2"},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                width=4 * scale,
                height=3 * scale
            )
        }
    )
    return g

def pseudotime_UMAP_dropdown():
    m = dcc.Dropdown(
        id='pseudotime_calculation_UMAP_dropdown',
        options=[
            {'label': 'leiden', 'value': 'leiden_n'},
            {'label': "# UMIs", "value": "total_counts"},
            {'label': "# UMIs [ln(1+UMIs)]", "value": "log1p_total_counts"},
            {'label': "# unique genes", "value": "n_genes"},
            {'label': 'psuedotime', 'value': 'pseudotime'},
            {'label': 'differentiation potential', 'value': 'differentiation potential'},
        ],
        value=None,
        placeholder="leiden",
        multi=False,
        searchable=True
        ) 
    return m 
