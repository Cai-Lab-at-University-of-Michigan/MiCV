import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import pandas as pd
import numpy as np
import anndata as ad

def plot_clustering_UMAP():
    g = dcc.Graph(
        id='Clustering_UMAP_plot',
        figure={
            'data': [
                go.Scattergl(
                    x=[0,1,2],
                    y=[-3,-2,-1],
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
                width=1024,
                height=800
            )
        }
    )
    return g

def old_clustering_UMAP(adata=None):
    g = dcc.Graph(
        id='Clustering_UMAP_plot',
        figure={
            'data': [
                go.Scattergl(
                    x=adata.obsm["X_umap"][adata.obs["leiden_n"] == i,0],
                    y=adata.obsm["X_umap"][adata.obs["leiden_n"] == i,1],
                    text=str(adata.obs["leiden_n"]),
                    mode='markers',
                    opacity=0.7,
                    marker={
                        "line_width": 1
                    },
                    name=str(i)
                ) for i in adata.obs["leiden_n"].unique()
            ],
            'layout': go.Layout(
                xaxis={'title': 'UMAP1'},
                yaxis={'title': "UMAP2"},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                width=1024,
                height=800
            )
        }
    )
    return g