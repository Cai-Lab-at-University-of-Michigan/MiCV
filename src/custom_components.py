'''
These functions generate placeholder figures/interactable components 
that get updated using callbacks after loading data.
'''

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import pandas as pd
import numpy as np
import anndata as ad

# scale of plot sizes
scale = 200


def plot_clustering_UMAP():
    g = dcc.Graph(
        id='clustering_UMAP_plot',
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

def plot_pseudotime_UMAP():
    g = dcc.Graph(
        id='pseudotime_UMAP_plot',
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
        id='pseudotime_gene_plot',
        figure={
            'data': [
                go.Scattergl(
                    x=None,
                    y=None,
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
                width=4 * scale,
                height=2 * scale
            )
        }
    )
    return g

def plot_expression_UMAP():
    g = dcc.Graph(
        id='expression_UMAP_plot',
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

def plot_gene_violin():
    g = dcc.Graph(
        id='violin_gene_plot',
        figure={
            'data': [
                go.Violin(
                    y=None,
                    text=str("NULL"),
                    opacity=0.7,
                    name=str("NULL"),
                    box_visible=True,
                    meanline_visible=True,
                    points="all"
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'Gene'},
                yaxis={'title': "Expression"},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                width=4 * scale,
                height=2 * scale
            )
        }
    )
    return g


def single_gene_dropdown():
    m = dcc.Dropdown(
        id='single_gene_dropdown',
        options=[
            {'label': 'GeneA', 'value': 'GeneA'}
        ],
        value=None,
        multi=False,
        searchable=True
        )  
    return m

def multi_gene_dropdown():
    m = dcc.Dropdown(
        id='multi_gene_dropdown',
        options=[
            {'label': 'GeneA', 'value': 'GeneA'}
        ],
        value=[],
        multi=True,
        searchable=True
        ) 
    return m 

def pseudotime_dropdown():
    m = dcc.Dropdown(
        id='pseudotime_dropdown',
        options=[
            {'label': 'pseudotime', 'value': 'pseudotime'},
            {'label': 'differentiation potential', 'value': 'differentiation_potential'},
        ],
        value=None,
        placeholder="pseudotime",
        multi=False,
        searchable=True
        ) 
    return m 

def clustering_dropdown():
    m = dcc.Dropdown(
        id='clustering_dropdown',
        options=[
            {'label': 'leiden', 'value': 'leiden_n'},
            {'label': 'user_0', 'value': 'user_0'},
            {'label': 'user_1', 'value': 'user_1'},
            {'label': 'user_2', 'value': 'user_2'},
            {'label': 'user_3', 'value': 'user_3'},
            {'label': 'user_4', 'value': 'user_4'},
        ],
        value=None,
        placeholder="leiden",
        multi=False,
        searchable=True
        ) 
    return m 

def gene_data_table():
    header = ["gene name", "flybase ID", 
              "human ortholog", "HGNC ID"]    
    data = ["null", "null", "null", "null"]
    
    t = html.Table(id="gene_data_table", children=[
        # Header
        html.Tr([html.Th(col) for col in header]),
        # Body
        html.Tr([html.Td(val) for val in data])
    ])
    return t

def plot_processing_UMAP():
    g = dcc.Graph(
        id='processing_UMAP_plot',
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

def processing_UMAP_dropdown():
    m = dcc.Dropdown(
        id='processing_UMAP_dropdown',
        options=[
            {'label': 'leiden', 'value': 'leiden_n'},
            {'label': 'psuedotime', 'value': 'pseudotime'},
            {'label': 'differentiation potential', 'value': 'differentiation potential'},
            {'label': "# UMIs (log1p)", "value": "# UMIs (log1p)"},
            {'label': "# unique genes", "value": "# unique genes"},
        ],
        value=None,
        placeholder="leiden",
        multi=False,
        searchable=True
        ) 
    return m 

def plot_processing_QC():
    g = dcc.Graph(
        id='processing_QC_plot',
        figure={
            'data': [
                go.Violin(
                    y=None,
                    text=str("NULL"),
                    opacity=0.7,
                    name=str("NULL"),
                    box_visible=True,
                    meanline_visible=True,
                    points="all"
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'QC factor'},
                yaxis={'title': "Counts"},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                width=4 * scale,
                height=2 * scale
            )
        }
    )
    return g

def processing_QC_dropdown():
    m = dcc.Dropdown(
    id='processing_QC_dropdown',
    options=[
        {'label': '# UMIs (log1p)', 'value': 'log1p_total_counts'},
        {'label': '# unique genes', 'value': 'n_genes'},
    ],
    value=None,
    placeholder="# UMIs (log1p)",
    multi=False,
    searchable=True
    ) 
    return m