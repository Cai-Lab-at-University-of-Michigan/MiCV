import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import plotly.graph_objs as go

import pandas as pd

from plotting.plotting_parameters import *

### PROCESSING TAB CUSTOM COMPONENTS ###
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
            ]
        }
    )
    return g

def processing_UMAP_dropdown():
    m = dcc.Dropdown(
        id='processing_UMAP_dropdown',
        options=[
            {'label': 'leiden', 'value': 'leiden_n'},
            {'label': "# UMIs ", "value": "total_counts"},
            {'label': "# UMIs [ln(1+UMIs)]", "value": "log1p_total_counts"},
            {'label': "# unique genes", "value": "n_genes"},
        ],
        value=None,
        placeholder="Select a cell observation",
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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                #width=4 * scale,
                #height=2 * scale
            )
        }
    )
    return g

def processing_QC_dropdown():
    m = dcc.Dropdown(
    id='processing_QC_dropdown',
    options=[
        {'label': "# UMIs ", "value": "total_counts"},
        {'label': "# UMIs [ln(1+UMIs)]", "value": "log1p_total_counts"},
        {'label': "# unique genes", "value": "n_genes"}
    ],
    value=None,
    placeholder="Select a QC parameter",
    multi=False,
    searchable=True
    ) 
    return m

def downsample_cells_slider():
    m = dbc.Card(children=[
            dbc.CardHeader("Downsample cells"),
            dbc.CardBody(children=[
                html.Div(children=[
                    html.P("Take only a fraction of cells for downstream analysis")
                ]),
                dcc.Slider(
                    id="downsample_cells_slider",
                    min=5,
                    max=100,
                    step=5,
                    marks={
                        10: "10%",
                        25: "25%",
                        50: "50%",
                        75: "75%",
                        100: "100% (default)"
                    },
                    value=100
                )
            ]),
            dbc.CardFooter(children=[
                html.H4([dbc.Badge("% Cells = 100% ()", 
                                   id="downsample_cells_output_container")
                ])
            ])
        ])
    return m

def downsample_counts_slider():
    m = dbc.Card(children=[
            dbc.CardHeader("Downsample UMI counts"),
            dbc.CardBody(children=[
                html.Div(children=[
                    html.P("Take only a fraction of UMIs (unique mRNA molecules) for downstream analysis")
                ]),
                dcc.Slider(
                    id="downsample_counts_slider",
                    min=5,
                    max=100,
                    step=5,
                    marks={
                        10: "10%",
                        25: "25%",
                        50: "50%",
                        75: "75%",
                        100: "100% (default)"
                    },
                    value=100
                )
            ]),
            dbc.CardFooter(children=[
                html.H4([dbc.Badge("% UMI counts = 100% ()", 
                                   id="downsample_counts_output_container")
                ])
            ])
        ])
    return m

def neighbors_method_radio():
    m = dbc.RadioItems(
        id="neighbors_method_radio",
        options=[
            {'label': 'standard', 'value': 'standard'},
            {'label': 'batch-corrected (bbknn)', 'value': 'bbknn'}
        ],
        value='standard'
        )
    return m

def n_dims_proj_radio():
    m = dbc.RadioItems(
        id="n_dims_proj_radio",
        options=[
            {'label': '2D ', 'value': 2},
            {'label': '3D ', 'value': 3}
        ],
        value=2
        )
    return m

def n_dims_processing_radio():
    m = dbc.RadioItems(
        id="n_dims_processing_radio",
        options=[
            {'label': '2D ', 'value': 2},
            {'label': '3D ', 'value': 3}
        ],
        value=2
        )
    return m
