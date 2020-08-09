import dash
from dash.dependencies import Input, Output, State

import os
import zipfile as zf
import base64
import io

from helper_functions import *
from plotting.plotting_functions import *
from status.status_functions import *
from app import app

from . processing_functions import *

@app.callback(
    [Output(f"{x}-collapse", "is_open") for x in ["downsample", "QC", "projection", "clustering"]],
    [Input(f"{x}-collapse-button", "n_clicks") for x in ["downsample", "QC", "projection", "clustering"]],
    [State(f"{x}-collapse", "is_open") for x in ["downsample", "QC", "projection", "clustering"]],
)
def toggle_procssing_accordion(n1, n2, n3, n4, is_open1, is_open2, is_open3, is_open4):
    ctx = dash.callback_context

    if not ctx.triggered:
        return dash.no_update
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if button_id == "downsample-collapse-button" and n1:
        return not is_open1, False, False, False
    elif button_id == "QC-collapse-button" and n2:
        return False, not is_open2, False, False
    elif button_id == "projection-collapse-button" and n3:
        return False, False, not is_open3, False
    elif button_id == "clustering-collapse-button" and n4:
        return False, False, False, not is_open4
    return False, False, False, False

@app.callback(
    Output('downsample_cells_output_container', 'children'),
    [Input('downsample_cells_slider', 'value')],
    [State('session-id', "children")]
)
def update_downsample_cells_slider_output(value, session_ID):
    state = cache_state(session_ID)
    if (state is None):
        return dash.no_update
    elif not ("# cells/obs" in state.keys()):
        return dash.no_update

    n_cells = state["# cells/obs"]
    return ("% Cells = " + str(value) + "% "
          + "(" + str(int((value/100) * n_cells)) + ")")

@app.callback(
    Output('downsample_counts_output_container', 'children'),
    [Input('downsample_counts_slider', 'value')],
    [State('session-id', "children")]
)
def update_downsample_counts_slider_output(value, session_ID):
    state = cache_state(session_ID)
    if (state is None):
        return dash.no_update
    elif not ("# counts" in state.keys()):
        return dash.no_update

    n_counts = state["# counts"]
    return ("% UMI counts = " + str(value) + "% "
          + "(" + str(int((value/100) * n_counts)) + ")")

@app.callback(
    Output('min_max_genes_slider_output_container', 'children'),
    [Input('min_max_genes_slider', 'value')]
)
def update_min_max_genes_output(value):
    return ("[min_genes, max_genes] = " + str(value))

@app.callback(
    Output('min_cells_slider_output_container', 'children'),
    [Input('min_cells_slider', 'value')]
)
def update_min_cells_output(value):
    return ("min_cells = " + str(value))

@app.callback(
    Output('n_top_genes_slider_output_container', 'children'),
    [Input('n_top_genes_slider', 'value')]
)
def update_n_top_genes_output(value):
    return ("# highly variable genes = " + str(value))

@app.callback(
    Output('n_neighbors_slider_output_container', 'children'),
    [Input('n_neighbors_slider', 'value')]
)
def update_n_neighbors_output(value):
    return ("n_neighbors = " + str(value))


@app.callback(
    Output('clustering_resolution_slider_output_container', 'children'),
    [Input('clustering_resolution_slider', 'value')]
)
def update_clustering_resolution_output(value):
    return "clustering resolution = " + str(value)

@app.callback(
    [Output("processing_UMAP_plot", "figure"),
     Output("refresh_all_status", "children")],
    [Input("refresh_all_button", "n_clicks"),
     Input("refresh_projection_button", "n_clicks"),
     Input("refresh_clustering_button", "n_clicks"),
     Input("processing_UMAP_dropdown", "value"),
     Input("neighbors_method_radio", "value"),
     Input("n_dims_processing_radio", "value")],
    [State("session-id", "children"),
     State("downsample_cells_slider", "value"),
     State("downsample_counts_slider", "value"),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value"),
     State("min_max_genes_slider", "value"),
     State("min_cells_slider", "value"),
     State("n_top_genes_slider", "value")]
)
def refresh_processing_UMAP(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, processing_plot_type,
                            neighborhood_method, n_dim_proj_plot, 
                            session_ID, pct_cells, pct_counts,
                            n_neighbors, resolution, 
                            min_max_genes, min_cells, n_top_genes,
                            adata=None, target_sum=1e6, 
                            flavor="cell_ranger", n_comps=50, random_state=0):
    
    default_return = [dash.no_update, dash.no_update]
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]


    if  (button_id == "refresh_clustering_button"):
        if (clust_btn_clicks in [None, 0]):
            return default_return
        
        n_steps = 2
        adata = cache_adata(session_ID)
        if (adata is None):
            return default_return
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        adata = do_clustering(session_ID, adata, resolution=resolution)
        cache_progress(session_ID, progress=int(2/n_steps * 100))

    elif(button_id == "refresh_projection_button"):
        if (proj_btn_clicks in [None, 0]):
            return default_return
        
        n_steps = 4
        adata = cache_adata(session_ID)
        if (adata is None):
            return default_return
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        adata, = do_PCA(session_ID, adata, n_comps=n_comps, 
                        random_state=random_state),   
        cache_progress(session_ID, progress=int(2/n_steps * 100))

        adata = do_neighborhood_graph(session_ID, adata, 
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state,
                                      method=neighborhood_method)
        cache_progress(session_ID, progress=int(3/n_steps * 100))
        
        adata = do_UMAP(session_ID, adata, random_state=random_state)
        cache_progress(session_ID, progress=int(4/n_steps * 100))

    elif(button_id == "refresh_all_button"):
        if (all_btn_clicks in [None, 0]):
            return default_return

        n_steps = 7
        adata = cache_adata(session_ID)
        if (adata is None):
            return default_return
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        adata = downsample_adata(session_ID, adata, pct_cells=pct_cells,
                                 pct_counts=pct_counts) 
        cache_progress(session_ID, progress=int(2/n_steps * 100))

        adata = preprocess_data(session_ID, adata, min_cells=min_cells,
                                min_genes=min_max_genes[0], 
                                max_genes=min_max_genes[1], 
                                target_sum=target_sum, flavor=flavor, 
                                n_top_genes=n_top_genes)
        cache_progress(session_ID, progress=int(3/n_steps * 100))

        adata, = do_PCA(session_ID, adata, n_comps=n_comps, 
                        random_state=random_state),
        cache_progress(session_ID, progress=int(4/n_steps * 100))

        adata = do_neighborhood_graph(session_ID, adata, neighborhood_method,
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state)
        cache_progress(session_ID, progress=int(5/n_steps * 100))

        adata = do_UMAP(session_ID, adata, random_state=random_state)
        cache_progress(session_ID, progress=int(6/n_steps * 100))

        adata = do_clustering(session_ID, adata, resolution=resolution)
        cache_progress(session_ID, progress=int(7/n_steps * 100))
        
        default_return[1] = "Processing successful"

    # if it's a dropdown menu update - load adata
    elif(button_id == "processing_UMAP_dropdown"
      or button_id == "n_dims_processing_radio"):
        if (processing_plot_type in [0, "", None, []]):
            return default_return
        else:
            if (adata_cache_exists(session_ID) is False):
                return default_return

    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return default_return
    
    # update the plot
    if (processing_plot_type in [0, "", None, []]):
        return default_return

    if (processing_plot_type == "leiden_n"):
        return plot_UMAP(session_ID, "leiden", n_dim=n_dim_proj_plot), "processing successful"
    elif (processing_plot_type == "pseudotime"):
        return plot_pseudotime_UMAP(session_ID, "pseudotime"), "processing successful"
    elif (processing_plot_type == "differentiation potential"):
        return plot_pseudotime_UMAP(session_ID, "differentiation_potential"), "processing successful"
    elif (processing_plot_type == "total_counts"):
        return plot_expression_UMAP(session_ID, "total_counts", n_dim=n_dim_proj_plot), "processing successful"
    elif (processing_plot_type == "n_genes"):
        return plot_expression_UMAP(session_ID, "n_genes", n_dim=n_dim_proj_plot), "processing successful"
    elif (processing_plot_type == "log1p_total_counts"):
        return plot_expression_UMAP(session_ID, "log1p_total_counts", n_dim=n_dim_proj_plot), "processing successful"

@app.callback(
    Output("processing_QC_plot", "figure"),
    [Input("processing_QC_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_violin_QC_plot(selected_QC, session_ID):
    default_return = dash.no_update
    if (selected_QC in [None, 0, []]):
        return default_return

    if (adata_cache_exists(session_ID) is False):
        return default_return

    # plot function expects list of factors/genes, but for QC
    # we will only show one at a time here - list is required though
    return plot_expression_violin(session_ID, [selected_QC], show_points = False)