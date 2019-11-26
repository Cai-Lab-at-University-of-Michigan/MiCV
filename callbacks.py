import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
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
     Input("clustering_dropdown", "value"),
     Input("Pseudotime_UMAP_plot", "selectedData"),
     Input("Expression_UMAP_plot", "selectedData"),
     Input("Violin_gene_plot", "selectedData")],
    [State("session-id", "children"),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value"),
     State("Clustering_UMAP_plot", "selectedData")]
)
def refresh_clustering_plot(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, load_btn_clicks, clustering_plot_type, pt_selected, expr_selected, violin_selected,
                            session_ID, n_neighbors, resolution, clust_selected,  
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

            if (clustering_plot_type in [0, "", None, []]):
                return dash.no_update



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
    
    elif(button_id in ["Pseudotime_UMAP_plot",
                       "Expression_UMAP_plot",
                       "Violin_gene_plot"]):
        if (pt_selected is None
        and expr_selected is None
        and violin_selected is None):
            return dash.no_update
        else:
            adata = cache_adata(session_ID, adata)



    # if it's a dropdown menu update - load adata
    elif(button_id == "clustering_dropdown"):
        if (clustering_plot_type in [0, "", None, []]):
            return dash.no_update
        else:
            adata = cache_adata(session_ID)

        # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return dash.no_update


    # on first run with this dataset, add columns for user cluster annotations
    for i in ["user_" + str(j) for j in range(0, 10)]:
        if not (i in adata.obs.columns):
            adata.obs[i] = ["unnassigned" for j in adata.obs.index.to_list()]

    # figure out which cells need to be selected, based on other graphs
    violin_selected = get_violin_intersection(session_ID, adata, violin_selected)
    selected_cell_intersection = get_cell_intersection(session_ID, adata,
                                                        [clust_selected, 
                                                         pt_selected,
                                                         expr_selected,
                                                         violin_selected])
    selected_points = []
    for c in selected_cell_intersection:
        selected_points.append((adata.obs).index.get_loc(c))
    if (len(selected_points) == 0):
        selected_points = range(0, len(((adata.obs).loc[:, "cell_ID"]).index.to_list()))
    
    # update the plot
    print("[STATUS] updating plot by: " + str(clustering_plot_type))
    traces = []
    for i in sorted(adata.obs[clustering_plot_type].unique()):
        a = adata[adata.obs[clustering_plot_type] == i]
        s = []
        for c in selected_cell_intersection:
            if (c in a.obs["cell_numeric_index"]):
                s.append((a.obs).index.get_loc(c))
        traces.append(
            go.Scattergl(
                x=a.obsm["X_umap"][:,0],
                y=a.obsm["X_umap"][:,1],
                text="Cell ID: " + a.obs["cell_ID"],
                mode='markers',
                selectedpoints=s,
                marker={
                    'size': 10,
                    'line': {'width': 1, 'color': 'grey'}
                },
                unselected={
                    "marker": {"opacity": 0.2,
                    }
                },
                selected={
                    "marker": {"opacity": 0.9,
                    }
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
        if (pt_btn_clicks in [None, 0, "", []]):
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
    traces.append(
        go.Scattergl(
            x=adata.obsm["X_umap"][:,0],
            y=adata.obsm["X_umap"][:,1],
            text="Cell ID: " + adata.obs["cell_ID"],
            mode='markers',
            marker={
                'size': 10,
                'line': {'width': 1, 'color': 'grey'},
                "color": adata.obs[str(pt_plot_type)],
                "colorscale": "plasma",
                "cmin": 0,
                "cmax": 1,
                "colorbar": dict(
                    title="Pseudotime"
                ),
            },
            unselected={
                "marker": {"opacity": 0.2,
                }
            },
            selected={
                "marker": {"opacity": 0.7,
                }
            },
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
    traces.append(
        go.Scattergl(
            x=adata.obsm["X_umap"][:,0],
            y=adata.obsm["X_umap"][:,1],
            text="Cell ID: " + adata.obs["cell_ID"],
            mode='markers',
            marker={
                'size': 10,
                'line': {'width': 1, 'color': 'grey'},
                "color": adata.raw.obs_vector(selected_gene),
                "colorscale": "viridis",
                "cmin": 0,
                "cmax": np.max(adata.raw.obs_vector(selected_gene)),
                "colorbar": dict(
                    title=str(selected_gene)
                ),
            },
            unselected={
                "marker": {"opacity": 0.2,
                }
            },
            selected={
                "marker": {"opacity": 0.7,
                }
            },
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
            transition = {'duration': 100},
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

    #adata = cache_adata(session_ID)
    gene_trends = cache_gene_trends(session_ID)

    # rearrange some data before plotting
    branches = list(gene_trends.keys())
    colors = pd.Series(sns.color_palette('Set2', len(selected_genes)).as_hex(), 
                       index=selected_genes)

    # regardless of what updates were requested - update the plot
    traces = []
    trends = gene_trends[branches[0]]["trends"]
    stds = gene_trends[branches[0]]["std"] * 25
    for i in selected_genes:
        if not (i in trends.index):
            print("[DEBUG] gene " + str(i)  + " not in gene trends; skipping")
            continue
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
    
    if (traces in [[], None]):
        print("[DEBUG] no traces added to violin plot")
        return dash.no_update

    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Pseudotime"},
            yaxis={"title": "Expression"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
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

    print("[STATUS] updating violin gene plot")
    adata = cache_adata(session_ID)

    # regardless of what updates were requested - update the plot
    traces = []
    for i in selected_genes:
        if not (i in adata.var.index):
            print("[DEBUG] gene " + str(i) + " not in var index; skipping")
            continue
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

    if (traces in [[], None]):
        print("[DEBUG] no traces added to violin plot")
        return dash.no_update

    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Gene"},
            yaxis={"title": "Expression"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
        )
    }

@app.callback(
    Output('define_cluster_text', 'children'),
    [Input('define_cluster_button', 'n_clicks')],
    [State('session-id', 'children'),
     State("clustering_dropdown", "value"),
     State("Clustering_UMAP_plot", "selectedData"),
     State("Pseudotime_UMAP_plot", "selectedData"),
     State("Expression_UMAP_plot", "selectedData"),
     State("Violin_gene_plot", "selectedData")]
)
def define_new_cluster(n_clicks, session_ID, clustering_group, 
                       clust_selected, pt_selected, expr_selected, violin_selected):
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
        i += 1

    # for violin plot selection, take the intersection of all points selected
    # i.e. if the user selected cells from multiple expression violin plots,
    # take the intersection of those cells (logical AND)
    violin_selected = get_violin_intersection(session_ID, adata, violin_selected) 
    
    # add the selected cells to this cluster
    selected_cells = get_cell_intersection(session_ID, adata, [clust_selected, 
                                                               pt_selected,
                                                               expr_selected,
                                                               violin_selected])
    
    (adata.obs).loc[selected_cells, clustering_group] = new_cluster_ID
    
    cache_adata(session_ID, adata)

    return ("put " + str(len(selected_cells)) + " cells in cluster " 
            + str(new_cluster_ID) + " under clustering group " + str(clustering_group))


@app.callback(
    Output('gene_data_table', 'children'),
    [Input('single_gene_dropdown', 'value')],
    [State('session-id', 'children')]
)
def update_gene_data_table(selected_gene, session_ID):
    print("[DEBUG] selected_gene: " + str(selected_gene))
    if (selected_gene in [None, "", 0, []]):
        return dash.no_update

    d = get_ortholog_data(session_ID, selected_gene)

    print("[DEBUG] d: " + str(d))
    
    snapshot = get_gene_snapshot(session_ID, selected_gene)
    # make the HTML table

    ret = []
    # gene snapshot
    snapshot_header = ["flybase ID", "gene symbol", "Flybase snapshot"]
    ret.append(html.Tr([html.Th(col) for col in snapshot_header]))
    items = [snapshot.iloc[0]["FBgn_ID"], 
             snapshot.iloc[0]["GeneSymbol"], 
             snapshot.iloc[0]["gene_snapshot_text"]]
    ret.append(html.Tr([html.Td(i) for i in items]))

    # human orthologs
    ortholog_header = ["human ortholog", "DIOPT score", "HGNC ID"]
    ret.append(html.Tr([html.Th(col) for col in ortholog_header]))
    if (len(d.index) == 0):
        ret.append(html.Tr([html.Td("N/A") for i in ortholog_header]))
    else:
        for index, row in d.iterrows():
            items = [row["Human_gene_symbol"],
                    row["DIOPT_score"],
                    row["Human_gene_HGNC_ID"]]
            items = [str(i) for i in items]
            r = html.Tr([html.Td(i) for i in items])
            ret.append(r)
    return ret