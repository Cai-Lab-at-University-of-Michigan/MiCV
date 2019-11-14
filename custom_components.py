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
                width=920,
                height=690
            )
        }
    )
    return g

def plot_pseudotime_UMAP():
    g = dcc.Graph(
        id='Pseudotime_UMAP_plot',
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
                width=920,
                height=690
            )
        }
    )
    return g

def gene_list_dropdown(gene_list):
    if (gene_list is None):
        gene_list = []
    return dcc.Dropdown(
        id='gene_list_dropdown',
        options=[{'label': i, 'value': i} for i in gene_list],
        value='Genes'
    )

def plot_gene_pseudotime():
    g = dcc.Graph(
        id='Pseudotime_gene_plot',
        figure={
            'data': [
                go.Scattergl(
                    x=[0,0.5,1],
                    y=[3,7,1],
                    text=str("NULL"),
                    opacity=0.7,
                    name=str("NULL"),
                    mode="lines+markers",
                    marker={
                        "opacity": 0,
                        "line_width": 2
                    }
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'Pseudotime'},
                yaxis={'title': "Expression"},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                width=920,
                height=300
            )
        }
    )
    return g

def plot_expression_UMAP():
    g = dcc.Graph(
        id='Expression_UMAP_plot',
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
                width=920,
                height=690
            )
        }
    )
    return g

def single_gene_dropdown():
    m = dcc.Dropdown(
        id='single_gene_dropdown',
        options=[
            {'label': 'GeneA', 'value': 'GeneA'},
            {'label': 'GeneB', 'value': 'GeneB'},
            {'label': 'GeneC', 'value': 'GeneC'}
        ],
        value=[],
        multi=False,
        searchable=True
        )  
    return m

def multi_gene_dropdown():
    m = dcc.Dropdown(
        id='multi_gene_dropdown',
        options=[
            {'label': 'GeneA', 'value': 'GeneA'},
            {'label': 'GeneB', 'value': 'GeneB'},
            {'label': 'GeneC', 'value': 'GeneC'}
        ],
        value=[],
        multi=True,
        searchable=True
        ) 
    return m 