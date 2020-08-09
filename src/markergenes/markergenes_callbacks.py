import pandas as pd

import dash
from dash.dependencies import Input, Output, State

from flask import send_file

import os
import io

from . markergenes_functions import identify_marker_genes

from plotting.plotting_functions import *
from helper_functions import *
from status.status_functions import *

from app import app

#### Marker gene page callbacks ####
@app.callback(
    Output("marker_gene_UMAP_dropdown", "options"),
    [Input("main_tabs", "active_tab")],
    [State('session-id', 'children')]
)
def refresh_marker_gene_UMAP_dropdown(active_tab, session_ID):
    default_return = dash.no_update

    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (active_tab != "markergenes_tab"):
        return default_return

    obs = cache_adata(session_ID, group="obs")
    if (obs is None):
        return default_return

    a = obs.select_dtypes(include=["category"])
    options = [
        {"label": str(x), "value": str(x)} for x in a.columns.to_list()
    ]
    return options

@app.callback(
    [Output("marker_gene_group_dropdown", "options"),
     Output("marker_gene_UMAP_plot", "figure")],
    [Input("marker_gene_UMAP_dropdown", "value"),
     Input("n_dims_proj_markergenes_radio", "value")],
    [State('session-id', 'children')]
)
def refresh_marker_gene_group_dropdown(obs_column,
                                       n_dims_proj, session_ID):
    default_return = [dash.no_update, dash.no_update]
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    obs = cache_adata(session_ID, group="obs")
    
    if (obs_column in ["", 0, None, []]):
        return default_return

    if ((obs is None)
    or not (obs_column in obs)):
        return default_return
    
    options = [
        {"label": str(x), "value": x} for x in sorted((obs[obs_column]).unique())
    ]
    options.insert(0, {"label": "all (default)", "value": "all"})
    
    if (obs_column in ["", None, 0, []]):
        fig = dash.no_update
    else:
        fig = plot_UMAP(session_ID, obs_column, n_dim=n_dims_proj)
    
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

    if (groups_to_rank in [[], None, 0, ""]):
        return dash.no_update
    
    n_steps = 3

    adata = cache_adata(session_ID)
    cache_progress(session_ID, progress=int(1/n_steps * 100))

    adata = identify_marker_genes(session_ID, adata, obs_column, groups_to_rank, method)
    cache_progress(session_ID, progress=int(2/n_steps * 100))

    cache_adata(session_ID, adata.uns, group="uns")
    generate_marker_gene_table(session_ID)
    cache_progress(session_ID, progress=int(3/n_steps * 100))

    return plot_marker_genes(session_ID, adata, obs_column, groups_to_rank)

@app.callback(
    [Output("download_marker_genes_button", "disabled"),
     Output("download_marker_genes_link", "href")],
    [Input("main_tabs", "active_tab"),
     Input("marker_gene_plot", "children")],
    [State("session-id", "children")]
)
def enabled_marker_genes_download_button(active_tab,
                                         plot_children, session_ID):
    default_return = [True, dash.no_update] #disabled
    
    if ((active_tab != "markergenes_tab")
    or  (plot_children in [0, None, "", []])):
        return default_return
    
    if (marker_genes_table_exists(session_ID)):
        return [False, "/download/" + session_ID + "/marker_genes"]
    else:
        return default_return

@app.server.route('/download/<path:path>/marker_genes')
def serve_marker_genes_table_csv(path):
    f = save_analysis_path + path + "/marker_genes.csv"
    if (os.path.isfile(f)):
        print("[DEBUG] file " + f + " found - serving")
        with open(f, "rb") as b:
            return send_file(io.BytesIO(b.read()), 
                             as_attachment=True,
                             attachment_filename="marker_genes.csv",
                             mimetype="application/octet-stream", cache_timeout=3
                    )
    else:
        print("[ERROR] file " + f + " not available for download")
        return ("Error - file not generated yet. Go back in your browser," 
             + "upload your raw data, and idenfity marker genes before downloading.")