from dash.dependencies import Input, Output, State
import plotly.graph_objs as go

from app import app
from analysis_functions import *

@app.callback(
    Output("Clustering_UMAP_plot", "figure"),
    [Input("refresh_all_button", "n_clicks")],
    [State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value")]
)
def refresh_everything(n_clicks, n_neighbors, resolution, data_dir=None,
                       min_cells=2, min_genes=200, max_genes=10000,
                       target_sum=1e6, flavor="cell_ranger", 
                       n_top_genes=2000, n_comps=50, random_state=0):
    
    print("[STATUS] refreshing everything")
    
    adata = preprocess_data(data_dir, min_cells=min_cells, min_genes=min_genes,
                            max_genes=max_genes, target_sum=target_sum,
                            flavor=flavor, n_top_genes=n_top_genes)
    
    adata, = do_PCA(adata, n_comps=n_comps, random_state=random_state),

    adata = do_neighborhood_graph(adata, n_neighbors=n_neighbors, 
                                  random_state=random_state)
    adata = do_UMAP(adata, random_state=random_state)
    
    adata = do_clustering(adata, resolution=resolution)
    
    print("[STATUS] updating plot")
    traces = []
    for i in adata.obs["leiden_n"].unique():
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