import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go

from plotting.plotting_parameters import *

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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                autosize=True
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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                autosize=True
            )
        }
    )
    return g

def plot_gene_pseudotime():
    g = dcc.Graph(
        id='pseudotime_gene_plot',
        figure={
            'data': [
                go.Scatter(
                    x=None,
                    y=None,
                    text=str("NULL"),
                    name=str("NULL"),
                    mode="markers+lines",
                    marker={
                        "opacity": 0,
                        "line_width": 2
                    }
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'Pseudotime'},
                yaxis={'title': "Expression"},
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                autosize=True
            )
        }
    )
    return g

def plot_expression_UMAP():
    PLOT_RESOLUTION:dict = dict(width=600,
                                height=500)
    g = dcc.Graph(
        config={"toImageButtonOptions": PLOT_RESOLUTION},
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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                autosize=True
            )
        }
    )
    return g


'''
                go.Violin(
                    y=None,
                    text=str("NULL"),
                    opacity=0.7,
                    name=str("NULL"),
                    box_visible=True,
                    meanline_visible=True,
                    points="all"
                )
'''
def plot_gene_violin():

    g = dcc.Graph(
        id='violin_gene_plot',
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
                xaxis={'title': 'Gene'},
                yaxis={'title': "Expression"},
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                autosize=True
            )
        }
    )
    return g

def single_gene_expression_radio():
    m = dbc.RadioItems(
        id="single_gene_expression_radio",
        options=[
            {'label': 'standard', 'value': 'standard'},
            {'label': 'mixed', 'value': 'mixed'}
        ],
        value='standard'
        )
    return m  

def single_gene_dropdown():
    m = dcc.Dropdown(
        id='single_gene_dropdown',
        options=[
        ],
        value=None,
        multi=False,
        searchable=True,
        placeholder="Select a gene (single, log expression)"
        )  
    return m

def mixed_gene_dropdown():
    m = dcc.Dropdown(
        id='mixed_gene_dropdown',
        options=[
        ],
        value=None,
        multi=True,
        searchable=True,
        placeholder="Select some genes (multiple, relative expression)"
        )  
    return m

def multi_gene_dropdown():
    m = dcc.Dropdown(
        id='multi_gene_dropdown',
        options=[
        ],
        value=[],
        multi=True,
        searchable=True,
        placeholder="Select some genes (multiple)"
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
        placeholder="(only if calculated previously) Choose a pseudotime observation",
        multi=False,
        searchable=True
        ) 
    return m 

def pseudotime_gene_relative_radio():
    m = dbc.RadioItems(
        id="pseudotime_gene_relative_radio",
        options=[
            {"label": "absolute", "value": "absolute"},
            {"label": "relative", "value": "relative"}
        ],
        value="absolute",
        )
    return m

def pseudotime_gene_branch_dropdown():
    m = dcc.Dropdown(
        id='pseudotime_gene_branch_dropdown',
        options=[
        ],
        value=None,
        multi=False,
        searchable=True,
        placeholder="(only if calculated previously) choose a pseudotime branch"
        ) 
    return m 


def clustering_dropdown():
    m = dcc.Dropdown(
        id='clustering_dropdown',
        options=[
            {'label': 'leiden', 'value': 'leiden'},
            {'label': 'user_0', 'value': 'user_0'},
            {'label': 'user_1', 'value': 'user_1'},
            {'label': 'user_2', 'value': 'user_2'},
            {'label': 'user_3', 'value': 'user_3'},
            {'label': 'user_4', 'value': 'user_4'},
        ],
        value=None,
        placeholder="Choose automatic (leiden) clustering or a user column to add your own clusters to",
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

def n_dims_proj_expression_radio():
    m = dbc.RadioItems(
        id="n_dims_proj_expression_radio",
        options=[
            {'label': '2D ', 'value': 2},
            {'label': '3D ', 'value': 3}
        ],
        value=2,
        inline=True
        )
    return m

def gene_violin_count():
    m = dbc.Badge(id="gene_violin_count", children=["# cells selected: 0"],
                  className="ml-1")
    return m

def gene_UMAP_count():
    m = dbc.Badge(id="gene_UMAP_count", children=["# cells selected: 0"],
                  className="ml-1")
    return m

def pseudotime_UMAP_count():
    m = dbc.Badge(id="pseudotime_UMAP_count", children=["# cells selected: 0"],
                  className="ml-1")
    return m

def clustering_UMAP_count():
    m = dbc.Badge(id="clustering_UMAP_count", children=["# cells selected: 0"],
                  className="ml-1")
    return m

def total_cell_count():
    m = html.Div(id="total_cell_count", children=[0],
                 style={"display": "none"})
    return m