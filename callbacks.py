import dash
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go

from app import app
from analysis_functions import *
from helper_functions import *
import custom_components as cc

@app.callback(
    Output("Clustering_UMAP_plot", "figure"),
    [Input("refresh_all_button", "n_clicks"),
     Input("refresh_projection_button", "n_clicks"),
     Input("refresh_clustering_button", "n_clicks")],
    [State('session-id', 'children'),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value")]
)
def refresh_clustering_plot(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, session_ID, n_neighbors, 
                            resolution, adata=None, data_dir=None, 
                            min_cells=2, min_genes=200, max_genes=10000,
                            target_sum=1e6, flavor="cell_ranger", 
                            n_top_genes=2000, n_comps=50, random_state=0):

    print("[STATUS] refreshing plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # get the adata object from the cache
    adata = cache_adata(session_ID)

    if  (button_id == "refresh_clustering_button"):
        if (clust_btn_clicks in [None, 0]):
            return dash.no_update
        adata = do_clustering(session_ID, adata, resolution=resolution)
    
    elif(button_id == "refresh_projection_button"):
        if (proj_btn_clicks in [None, 0]):
            return dash.no_update   
        adata = do_neighborhood_graph(session_ID, adata, 
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state)
        adata = do_UMAP(session_ID, adata, random_state=random_state)
    
    elif(button_id == "refresh_all_button"):
        if (all_btn_clicks in [None, 0]):
            return dash.no_update
        print("[STATUS] refreshing everything")   
        adata = preprocess_data(session_ID, data_dir, min_cells=min_cells,
                                min_genes=min_genes, max_genes=max_genes, 
                                target_sum=target_sum, flavor=flavor, 
                                n_top_genes=n_top_genes)
        adata, = do_PCA(session_ID, adata, n_comps=n_comps, 
                        random_state=random_state),
        adata = do_neighborhood_graph(session_ID, adata, 
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state)
        adata = do_UMAP(session_ID, adata, random_state=random_state)
        adata = do_clustering(session_ID, adata, resolution=resolution)
    
    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return dash.no_update
    
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating plot")
    traces = []
    for i in sorted(adata.obs["leiden_n"].unique()):
        a = adata[adata.obs["leiden_n"] == i]
        traces.append(
            go.Scattergl(
                x=a.obsm["X_umap"][:,0],
                y=a.obsm["X_umap"][:,1],
                text="Cell ID: " + a.obs["cell_ID"],
                mode='markers',
                opacity=0.7,
                marker={
                    'size': 10,
                    'line': {'width': 1, 'color': 'grey'}
                },
                name=("Cluster " + str(i))
            )
        )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "UMAP 1"},
            yaxis={"title": "UMAP 2"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
        )
    }

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
    Output("Pseudotime_UMAP_plot", "figure"),
    [Input("refresh_pseudotime_button", "n_clicks")],
    [State('session-id', 'children'),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value")]
)
def refresh_pseudotime_plot(pt_btn_clicks, session_ID, n_neighbors, 
                            resolution, adata=None, data_dir=None):

    print("[STATUS] refreshing plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # get the adata object from the cache
    adata = cache_adata(session_ID)

    if  (button_id == "refresh_pseudotime_button"):
        if (pt_btn_clicks in [None, 0]):
            return dash.no_update
        adata = do_pseudotime(session_ID, adata)
    
    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return dash.no_update
    
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating pseudotime plot")
    traces = []
    for i in sorted(adata.obs["leiden_n"].unique()):
        a = adata[adata.obs["leiden_n"] == i]
        traces.append(
            go.Scattergl(
                x=a.obsm["X_umap"][:,0],
                y=a.obsm["X_umap"][:,1],
                text="Cell ID: " + a.obs["cell_ID"],
                mode='markers',
                opacity=0.7,
                marker={
                    'size': 10,
                    'line': {'width': 1, 'color': 'grey'},
                    "color": a.obs["pseudotime"],
                    "colorscale": "plasma"
                },
                name=("Cluster " + str(i))
            )
        )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "UMAP 1"},
            yaxis={"title": "UMAP 2"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
        )
    }