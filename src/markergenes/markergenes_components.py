import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go

from plotting.plotting_parameters import *

### MARKER GENE TAB CUSTOM COMPONENTS ###
def marker_gene_UMAP_dropdown():
    m = dcc.Dropdown(
    id='marker_gene_UMAP_dropdown',
    options=[
    ],
    value=None,
    placeholder="Choose a cellular observation",
    multi=False,
    searchable=True
    ) 
    return m

def plot_marker_gene_UMAP():
    g = dcc.Graph(
        id='marker_gene_UMAP_plot',
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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                #width=4 * scale,
                #height=3 * scale
            )
        }
    )
    return g

def marker_gene_group_dropdown():
    m = dcc.Dropdown(
        id="marker_gene_group_dropdown",
        options=[
        ],
        value=None,
        placeholder="Choose some groups to compare",
        multi=True,
        searchable=True
        ) 
    return m 

def marker_gene_method_radio():
    m = dbc.RadioItems(
        id="marker_gene_method_radio",
        options=[
            {'label': 'logreg', 'value': 'logreg'},
            {'label': 't-test', 'value': 't-test'},
            {'label': 't-test_overestim_var', 'value': 't-test_overestim_var'},
            {'label': 'wilcoxon', 'value': 'wilcoxon'},
        ],
        value='logreg'
        )
    return m

def marker_gene_plot():
    m = html.Div(
        id="marker_gene_plot",
        children=[
            html.Img(src=None)
        ])
    return m

def n_dims_proj_markergenes_radio():
    m = dbc.RadioItems(
        id="n_dims_proj_markergenes_radio",
        options=[
            {'label': '2D ', 'value': 2},
            {'label': '3D ', 'value': 3}
        ],
        value=2
        )
    return m