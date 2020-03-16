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

def processing_dataset_dropdown():
    m = dbc.Col(children=[
            dbc.Row(children=[
                dbc.Col(
                    dcc.Dropdown(
                        id='processing_dataset_dropdown',
                        options=[
                            {'label': 'Cocanougher et. al. (2020) whole fly CNS', 'value': "00002"},
                            {'label': 'Davie et. al. (2018) aging fly brain', 'value': "00003"},
                            {'label': '10Xv3 5K PBMC', 'value': "00004"}, 
                            {'label': 'Sharma et. al. (2020) mouse somatosensory neurons', 'value': "00005"},
                            {'label': 'Zeisel et. al. (2018) mouse nervous system', 'value': "00006"},
                        ],
                        value=None,
                        placeholder="Select a pre-made dataset",
                        searchable=True,
                        multi=False,
                    ),
                width=9),
                dbc.Col(children=[
                    html.Div(children=[
                        dbc.Button(
                            children=["Load selected dataset"],
                            id="processing_load_dataset_button",
                        ),
                        dcc.Loading(children=[
                            html.Div(children="", 
                                     id='load_selected_dataset_success_output',
                                     style={'margin-top': 10}
                            )
                        ])
                    ])
                ], width=3),
            ]),
        ], width=8)
    return m

def processing_data_upload():
    m = dbc.Col(children=[
            dcc.Upload(
                id='upload_raw_data',
                children=html.Div([
                    '''Drag and drop either an h5ad anndata object,
                    or a .zip file containing your 10X output
                    directory's contents. Alternatively, ''',
                    html.A('click here to select a file.')
                ]),
                style={
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
            html.Div(id='upload_raw_data_success_output',
                     style={'margin-top': 20})
        ], width=4)
    return m