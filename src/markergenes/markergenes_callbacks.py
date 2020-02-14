import pandas as pd

import dash
from dash.dependencies import Input, Output, State

from . markergenes_functions import identify_marker_genes

from plotting_functions import *
from helper_functions import *
from app import app

#### Marker gene page callbacks ####
@app.callback(
    Output("marker_gene_UMAP_dropdown", "options"),
    [Input("refresh_all_status", "children")],
    [State('session-id', 'children')]
)
def refresh_marker_gene_UMAP_dropdown(processing_status, session_ID):
    print("[STATUS] refreshing marker gene UMAP dropdown options")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if  (True):
        print("[DEBUG] processing_status: " + str(processing_status))
        if (processing_status in [None, [], 0, ""]):
            return dash.no_update

    adata = cache_adata(session_ID)
    if (adata is None):
        return dash.no_update

    a = adata.obs.select_dtypes(include=["category"])
    options = [
        {"label": str(x), "value": str(x)} for x in a.columns.to_list()
    ]
    return options

@app.callback(
    [Output("marker_gene_group_dropdown", "options"),
     Output("marker_gene_UMAP_plot", "figure")],
    [Input("marker_gene_UMAP_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_marker_gene_group_dropdown(obs_column, session_ID):
    print("[STATUS] refreshing marker gene group dropdown options")
    default_return = [dash.no_update, dash.no_update]
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return [dash.no_update, dash.no_update]
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if not (button_id == "marker_gene_UMAP_dropdown"): 
        return default_return

    print("[DEBUG] loading adata")
    adata = cache_adata(session_ID)
    
    if ((adata is None) or (obs_column is None)):
        return default_return
    options = [
        {"label": str(x), "value": x} for x in (adata.obs[obs_column]).unique()
    ]
    options.insert(0, {"label": "all (default)", "value": "all"})
    fig = plot_UMAP(adata, obs_column)
    
    return [options, fig]

@app.callback(
    Output("marker_gene_plot", "children"),
    [Input("recalculate_marker_genes", "n_clicks")],
    [State("marker_gene_UMAP_dropdown", "value"),
     State("marker_gene_group_dropdown", "value"),
     State("marker_gene_method_radio", "value"),
     State('session-id', 'children')]
)
def refresh_marker_gene_plot(n_clicks, obs_column, groups_to_rank, 
                             method, session_ID):
    print("[STATUS] refreshing marker gene plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if not (button_id == "recalculate_marker_genes"):
        return dash.no_update

    if (n_clicks in [[], None, 0]):
        return dash.no_update

    adata = cache_adata(session_ID)
    adata = identify_marker_genes(adata, obs_column, groups_to_rank, method)
    return plot_marker_genes(adata, obs_column, groups_to_rank)