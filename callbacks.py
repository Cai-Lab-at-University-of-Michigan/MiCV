import dash
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import time
import pandas as pd
import seaborn as sns

from app import app
from analysis_functions import *
from helper_functions import *
import custom_components as cc

@app.callback(
    Output("Clustering_UMAP_plot", "figure"),
    [Input("refresh_all_button", "n_clicks"),
     Input("refresh_projection_button", "n_clicks"),
     Input("refresh_clustering_button", "n_clicks"),
     Input("load_analysis_button", "n_clicks")],
    [State('session-id', 'children'),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value")]
)
def refresh_clustering_plot(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, load_btn_clicks,
                            session_ID, n_neighbors, resolution, adata=None, 
                            data_dir=None, min_cells=2, min_genes=200, 
                            max_genes=10000, target_sum=1e6, flavor="cell_ranger", 
                            n_top_genes=2000, n_comps=50, random_state=0):

    print("[STATUS] refreshing plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Must go first - otherwise there's no adata object to load from cache
    if(button_id == "load_analysis_button"):
        if (load_btn_clicks in [None, 0]):
            return dash.no_update
        else:
            adata = sc.read_h5ad(save_analysis_path + "adata_cache.h5ad")
            #print(adata)
            cache_adata(session_ID, adata)

            adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])
            adata.obs["cell_ID"] = adata.obs.index
            adata.obs["cell_numeric_index"] = [i for i in range(0,len(adata.obs.index))]
            gene_trends = read_gene_trends(session_ID)
            print(gene_trends)
            #adata.uns["gene_trends"] = gene_trends
            cache_gene_trends(session_ID, gene_trends)  



    # get the adata object from the cache
    #adata = cache_adata(session_ID)

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
        #print(adata.uns.keys())
        if not ("pseudotime" in adata.obs):
            adata = do_pseudotime(session_ID, adata)
        #adata = do_pseudotime(session_ID, adata)
    
    
    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return dash.no_update
    
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating pseudotime plot")
    traces = []
    n = 0
    for i in sorted(adata.obs["leiden_n"].unique()):
        n += 1
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
                    "colorscale": "plasma",
                    "cmin": 0,
                    "cmax": 1,
                    "colorbar": dict(
                        title="Pseudotime"
                    ) if (n == len(sorted(adata.obs["leiden_n"].unique()))) 
                      else None,
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



save_analysis_path = "/home/nigelmic/data/programming/dash/scanpy/cache/"
@app.callback(
    Output("null_container_0", "children"),
    [Input("save_analysis_button", "n_clicks")],
    [State('session-id', 'children')]
)
def save_analysis(save_btn_clicks, session_ID):
    
    # figure out which button was triggered
    ctx = dash.callback_context
    if (not ctx.triggered):
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            
    if(button_id == "save_analysis_button"):
        if (save_btn_clicks in [None, 0]):
            pass
        else:
            adata = cache_adata(session_ID)
            adata.write(save_analysis_path + "adata_cache.h5ad")
        return dash.no_update

@app.callback(
    Output("single_gene_dropdown", "options"),
    [Input("refresh_all_button", "n_clicks"),
     Input("load_analysis_button", "n_clicks")],
    [State("session-id", "children")]
)
def update_single_gene_dropdown(n0, n1, session_ID):
    adata = cache_adata(session_ID)
    if (adata is None):
        i = 0
        while(i < 10):
            adata = cache_adata(session_ID)
            time.sleep(1)
            i += 1
        if (adata is None):
            return dash.no_update
    all_genes = adata.var.index.tolist()
    all_genes.sort()
    options = [{"label": i, "value": i} for i in all_genes]
    return options


@app.callback(
    Output("Expression_UMAP_plot", "figure"),
    [Input("single_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_expression_UMAP_plot(selected_gene, session_ID):
    if (selected_gene in [None, 0, []]):
        return dash.no_update

    print(selected_gene)
    adata = cache_adata(session_ID)
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating expression UMAP plot")
    traces = []
    n = 0 # counts through the list of leiden_n
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
                    "color": a.raw.obs_vector(selected_gene),
                    "colorscale": "viridis",
                    "cmin": 0,
                    "cmax": np.max(adata.raw.obs_vector(selected_gene)),
                    "colorbar": dict(
                        title=str(selected_gene)
                    ) if (n == len(sorted(adata.obs["leiden_n"].unique())) - 1) 
                      else None,
                },
                name=("Cluster " + str(i))
            )
        )
        n += 1
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
    Output("multi_gene_dropdown", "options"),
    [Input("refresh_all_button", "n_clicks"),
     Input("load_analysis_button", "n_clicks")],
    [State("session-id", "children")]
)
def update_multi_gene_dropdown(n0, n1, session_ID):
    adata = cache_adata(session_ID)
    if (adata is None):
        i = 0
        while(i < 10):
            adata = cache_adata(session_ID)
            time.sleep(1)
            i += 1
        if (adata is None):
            return dash.no_update
    all_genes = adata.var.index.tolist()
    all_genes.sort()
    options = [{"label": i, "value": i} for i in all_genes]
    return options

@app.callback(
    Output("Pseudotime_gene_plot", "figure"),
    [Input("multi_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_pseudotime_gene_plot(selected_genes, session_ID):
    if (selected_genes in [None, 0, []]):
        return dash.no_update

    print(selected_genes)
    adata = cache_adata(session_ID)
    gene_trends = read_gene_trends(session_ID)
    #gene_trends = adata.uns["gene_trends"]
    print(gene_trends)

    # rearrange some data before plotting
    branches = list(gene_trends.keys())
    colors = pd.Series(sns.color_palette('Set2', len(selected_genes)).as_hex(), 
                       index=selected_genes)
    print("[DEBUG] branches: " + str(branches) + "\n gene_trends: " + str(gene_trends))

    # regardless of what updates were requested - update the plot
    print("[STATUS] updating pseudotime gene trend plot")
    print("[DEBUG] gene_trends " + str(gene_trends))
    traces = []
    for i in selected_genes:
        trends = gene_trends[branches[0]]["trends"]
        stds = gene_trends[branches[0]]["std"]
        traces.append(
            go.Scattergl(
                x=trends.columns,
                y=trends.loc[i, :],
                text=str(i),
                mode="lines+markers",
                opacity=0.7,
                marker={
                    'size': 5,
                    'line': {'width': 2, 'color': colors[i]},
                    "color": colors[i],
                    "opacity": 0.25
                },
                name=(str(i))
            )
        )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Pseudotime"},
            yaxis={"title": "Expression"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
        )
    }



@app.callback(
    Output("Violin_gene_plot", "figure"),
    [Input("multi_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_violin_gene_plot(selected_genes, session_ID):
    if (selected_genes in [None, 0, []]):
        return dash.no_update

    print(selected_genes)
    adata = cache_adata(session_ID)

    # regardless of what updates were requested - update the plot
    print("[STATUS] updating violin gene plot")
    traces = []
    for i in selected_genes:
        d = adata.raw.obs_vector(i)
        traces.append(
            go.Violin(
                y=d,
                text="Cell ID: " + adata.obs["cell_ID"],
                opacity=0.7,
                box_visible=True,
                meanline_visible=True,
                points="all",
                name=(str(i))
            )
        )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Gene"},
            yaxis={"title": "Expression"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
        )
    }