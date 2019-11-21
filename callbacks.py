import dash
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import time
import re #regex
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
     Input("load_analysis_button", "n_clicks"),
     Input("clustering_dropdown", "value")],
    [State('session-id', 'children'),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value")]
)
def refresh_clustering_plot(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, load_btn_clicks, clustering_plot_type,
                            session_ID, n_neighbors, resolution, 
                            adata=None, data_dir=None, min_cells=2, min_genes=200, 
                            max_genes=10000, target_sum=1e6, flavor="cell_ranger", 
                            n_top_genes=2000, n_comps=50, random_state=0):

    print("[STATUS] refreshing clustering plot")
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
            adata = cache_adata(session_ID)

            adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])
            adata.obs["cell_ID"] = adata.obs.index
            adata.obs["cell_numeric_index"] = [i for i in range(0,len(adata.obs.index))]
            gene_trends = cache_gene_trends(session_ID)  

            gene_list = adata.var.index.tolist()
            cache_gene_list(session_ID, gene_list)



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

    # if it's a dropdown menu update - load adata
    else:
        adata = cache_adata(session_ID)

    
    if (clustering_plot_type in [0, "", None, []]):
        return dash.no_update

    for i in ["user_" + str(j) for j in range(0, 10)]:
        if not (i in adata.obs.columns):
            adata.obs[i] = ["unnassigned" for j in adata.obs.index.to_list()]
    adata = cache_adata(session_ID, adata)

    # update the plot
    print("[STATUS] updating plot by: " + str(clustering_plot_type))
    traces = []
    for i in sorted(adata.obs[clustering_plot_type].unique()):
        a = adata[adata.obs[clustering_plot_type] == i]
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
    [Input("refresh_pseudotime_button", "n_clicks"),
     Input("pseudotime_dropdown", "value")],
    [State('session-id', 'children'),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value")]
)
def refresh_pseudotime_plot(pt_btn_clicks, pt_plot_type, session_ID, n_neighbors, 
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
        if not (("pseudotime" in adata.obs)
            and ("differentiation_potential" in adata.obs)):
            adata = do_pseudotime(session_ID, adata)
        #adata = do_pseudotime(session_ID, adata)
    
    
    # do nothing if no buttons pressed
    elif((button_id == "not_triggered")
      or (pt_plot_type in [0, "", None, []])):
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
                    "color": a.obs[str(pt_plot_type)],
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
            cache_adata(session_ID, adata)
        return dash.no_update

@app.callback(
    Output("single_gene_dropdown", "options"),
    [Input("refresh_all_button", "n_clicks"),
     Input("load_analysis_button", "n_clicks")],
    [State("session-id", "children")]
)
def update_single_gene_dropdown(n0, n1, session_ID):
    if (n0 in [0, None] and n1 in [0, None]):
        return dash.no_update
    gene_list = cache_gene_list(session_ID)
    gene_list.sort()
    options = [{"label": i, "value": i} for i in gene_list]
    return options


@app.callback(
    Output("Expression_UMAP_plot", "figure"),
    [Input("single_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_expression_UMAP_plot(selected_gene, session_ID):
    if (selected_gene in [None, 0, []]):
        return dash.no_update

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
    if (n0 in [0, None] and n1 in [0, None]):
        return dash.no_update
    gene_list = cache_gene_list(session_ID)
    gene_list.sort()
    options = [{"label": i, "value": i} for i in gene_list]
    return options

@app.callback(
    Output("Pseudotime_gene_plot", "figure"),
    [Input("multi_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_pseudotime_gene_plot(selected_genes, session_ID):
    if (selected_genes in [None, 0, []]):
        return dash.no_update

    adata = cache_adata(session_ID)
    gene_trends = cache_gene_trends(session_ID)

    # rearrange some data before plotting
    branches = list(gene_trends.keys())
    colors = pd.Series(sns.color_palette('Set2', len(selected_genes)).as_hex(), 
                       index=selected_genes)

    # regardless of what updates were requested - update the plot
    print("[STATUS] updating pseudotime gene trend plot")
    traces = []
    for i in selected_genes:
        trends = gene_trends[branches[0]]["trends"]
        stds = gene_trends[branches[0]]["std"] * 25
        traces.append(
            go.Scattergl(
                x=trends.columns,
                y=trends.loc[i, :],
                text=str(i),
                mode="lines+markers",
                opacity=0.7,
                marker={
                    'size': stds.loc[i, :],
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

@app.callback(
    Output('define_cluster_text', 'children'),
    [Input('define_cluster_button', 'n_clicks')],
    [State('session-id', 'children'),
     State("clustering_dropdown", "value"),
     State("Clustering_UMAP_plot", "selectedData")]
)
def define_new_cluster(n_clicks, session_ID, clustering_group, selectedData):
    if (n_clicks in [0, None, "", []]):
        return dash.no_update

    if (clustering_group == "leiden_n"):
        return "will not add cells to automatically generated group [leiden]"

    adata = cache_adata(session_ID)
    adata.obs[clustering_group] = adata.obs[clustering_group].astype(str)

    # create a new unique cluster ID for this cluster in this clustering_group
    previous_cluster_IDs = adata.obs[clustering_group].unique()
    new_cluster_ID = str(len(previous_cluster_IDs))
    i = 0
    while (new_cluster_ID in previous_cluster_IDs):
        new_cluster_ID = str(len(previous_cluster_IDs) + i) 

    # add the selected cells to this cluster
    selected_cells = []
    for cell in selectedData["points"]:
        cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
        selected_cells.append(cell_ID)
        print("[DEBUG] cell_ID: " + cell_ID)
        #else:
        #    print("[DEBUG] no match for: " + cell["text"])

    
    (adata.obs).loc[selected_cells, clustering_group] = new_cluster_ID
    
    cache_adata(session_ID, adata)

    return ("put " + str(len(selected_cells)) + " cells in cluster " 
            + str(new_cluster_ID) + " under clustering group " + str(clustering_group))