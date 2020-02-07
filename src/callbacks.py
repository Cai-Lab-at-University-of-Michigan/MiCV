import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import plotly.graph_objs as go
import time
import re #regex
import pandas as pd
import seaborn as sns

import os
import zipfile as zf
import base64
import io

from app import app
from analysis_functions import *
import helper_functions
from helper_functions import *
import custom_components as cc
import plotting_functions as plotting

#### Analysis page callbacks ####

@app.callback(
    Output("clustering_UMAP_plot", "figure"),
    [Input("load_analysis_button", "n_clicks"),
     Input("define_cluster_button", "n_clicks"),
     Input("clustering_dropdown", "value"),
     Input("pseudotime_UMAP_plot", "selectedData"),
     Input("expression_UMAP_plot", "selectedData"),
     Input("violin_gene_plot", "selectedData"),
     Input("pseudotime_gene_plot", "selectedData")],
    [State("session-id", "children"),
     State("clustering_UMAP_plot", "selectedData")]
)
def refresh_clustering_plot(load_btn_clicks, 
                            cluster_btn_clicks, clustering_plot_type, pt_selected, expr_selected, violin_selected, pt_gene_selected, 
                            session_ID, clust_selected,  
                            adata=None, data_dir=None):

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
            adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
            cache_adata(session_ID)

            gene_trends = cache_gene_trends(session_ID)  


            gene_list = adata.var.index.tolist()
            cache_gene_list(session_ID, gene_list)

            if (clustering_plot_type in [0, "", None, []]):
                return dash.no_update
    
    elif(button_id == "define_cluster_button"):
        if (cluster_btn_clicks in [0, None, "", []]):
            return dash.no_update

        adata = cache_adata(session_ID)
        clustering_group = clustering_plot_type
        if not (clustering_group == "leiden_n"):
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
            
            pt_min, pt_max = get_pseudotime_min_max(session_ID, pt_gene_selected)

            # add the selected cells to this cluster
            selected_cells = get_cell_intersection(session_ID, 
                             adata, [clust_selected, pt_selected,
                                     expr_selected, violin_selected],
                             pt_min, pt_max)
    
            (adata.obs).loc[selected_cells, clustering_group] = new_cluster_ID
            
            cache_adata(session_ID, adata)

    elif(button_id in ["Pseudotime_UMAP_plot",
                       "Expression_UMAP_plot",
                       "Violin_gene_plot",
                       "Pseudotime_gene_plot"]):
        if (pt_selected is None
        and expr_selected is None
        and violin_selected is None
        and pt_gene_selected is None):
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
    pt_min, pt_max = get_pseudotime_min_max(session_ID, pt_gene_selected)
    selected_cell_intersection = get_cell_intersection(session_ID, adata,
                                                        [clust_selected, 
                                                         pt_selected,
                                                         expr_selected,
                                                         violin_selected],
                                                         pt_min, pt_max)
    selected_points = []
    for c in selected_cell_intersection:
        selected_points.append((adata.obs).index.get_loc(c))
    if (len(selected_points) == 0):
        selected_points = range(0, len(((adata.obs).loc[:, "cell_ID"]).index.to_list()))
    
    # update the plot
    return plotting.plot_UMAP(adata, clustering_plot_type, selected_cell_intersection)


@app.callback(
    Output("pseudotime_UMAP_plot", "figure"),
    [Input("pseudotime_dropdown", "value")],
    [State('session-id', 'children'),
     State("load_analysis_button", "n_clicks")]
)
def refresh_pseudotime_plot(pt_plot_type, session_ID, load_btn_clicks,
                            adata=None, data_dir=None):

    print("[STATUS] refreshing pseudotime plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (button_id == "pseudotime_dropdown"):
        if (load_btn_clicks in [None, "", [], 0]):
            return dash.no_update

    # do nothing if no buttons pressed
    elif((button_id == "not_triggered")
      or (pt_plot_type in [0, "", None, []])):
        return dash.no_update
    
    # get the adata object from the cache
    adata = cache_adata(session_ID)
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating pseudotime plot")
    return plotting.plot_pseudotime_UMAP(adata, pt_plot_type)



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
    [Input("load_analysis_button", "n_clicks")],
    [State("session-id", "children")]
)
def update_single_gene_dropdown(n0, session_ID):
    if (n0 in [0, None]):
        return dash.no_update
    gene_list = cache_gene_list(session_ID)
    gene_list.sort()
    options = [{"label": i, "value": i} for i in gene_list]
    return options


@app.callback(
    Output("expression_UMAP_plot", "figure"),
    [Input("single_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_expression_UMAP_plot(selected_gene, session_ID):
    if (selected_gene in [None, 0, []]):
        return dash.no_update

    adata = cache_adata(session_ID)
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating expression UMAP plot")
    return plotting.plot_expression_UMAP(adata, selected_gene)


@app.callback(
    Output("multi_gene_dropdown", "options"),
    [Input("load_analysis_button", "n_clicks")],
    [State("session-id", "children")]
)
def update_multi_gene_dropdown(n0, session_ID):
    if (n0 in [0, None]):
        return dash.no_update
    gene_list = cache_gene_list(session_ID)
    gene_list.sort()
    options = [{"label": i, "value": i} for i in gene_list]
    return options

@app.callback(
    Output("pseudotime_gene_plot", "figure"),
    [Input("multi_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_pseudotime_gene_plot(selected_genes, session_ID):
    print("[STATUS] plotting gene pseudotime expression for: " + str(selected_genes))
    if (selected_genes in [None, 0, []]):
        return dash.no_update

    #adata = cache_adata(session_ID)
    gene_trends = cache_gene_trends(session_ID)

    # rearrange some data before plotting
    branches = list(gene_trends.keys())

    # regardless of what updates were requested - update the plot
    # TODO: intelligent branch selection
    return plotting.plot_expression_trend(gene_trends, branches[1], selected_genes)



@app.callback(
    Output("violin_gene_plot", "figure"),
    [Input("multi_gene_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_violin_gene_plot(selected_genes, session_ID):
    if (selected_genes in [None, 0, []]):
        return dash.no_update

    print("[STATUS] updating violin gene plot")
    adata = cache_adata(session_ID)

    # regardless of what updates were requested - update the plot
    return plotting.plot_expression_violin(adata, selected_genes)

@app.callback(
    Output('gene_data_table', 'children'),
    [Input('single_gene_dropdown', 'value')],
    [State('session-id', 'children')]
)
def update_gene_data_table(selected_gene, session_ID):
    print("[STATUS] getting data for selected_gene: " 
          + str(selected_gene))

    if (selected_gene in [None, "", 0, []]):
        return dash.no_update

    d = get_ortholog_data(session_ID, selected_gene)
    d = d.sort_values(by="DIOPT_score", ascending=False, na_position="last")

    snapshot = get_gene_snapshot(session_ID, selected_gene)

    disease = get_disease_data(session_ID, selected_gene)
    
    # make the HTML table

    ret = []
    # gene snapshot
    snapshot_header = ["flybase ID", "gene symbol", "Flybase snapshot"]
    ret.append(html.Tr([html.Th(col) for col in snapshot_header]))
    if (len(snapshot.index) == 0):
        ret.append(html.Tr([html.Td("N/A") for i in snapshot_header]))
    else:
        items = [snapshot.iloc[0]["FBgn_ID"], 
                 snapshot.iloc[0]["GeneSymbol"], 
                 snapshot.iloc[0]["gene_snapshot_text"]]
        ret.append(html.Tr([html.Td(html.A(items[0], 
                                           href="http://flybase.org/reports/" 
                                           + items[0] + ".html",
                                           target="_blank")),
                            html.Td(items[1]),
                            html.Td(items[2])
                            ]))

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
            r = html.Tr([html.Td(items[0]),
                         html.Td(items[1]),
                         html.Td(html.A(items[2], 
                                        href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"
                                        + items[2],
                                        target="_blank"))
                        ])
            ret.append(r)

    # disease data
    disease_header = ["human disease context", "based on orthology with",
                      "flybase reference"]
    ret.append(html.Tr([html.Th(col) for col in disease_header]))
    if (len(disease.index) == 0):
        ret.append(html.Tr([html.Td("N/A") for i in disease_header]))
    else:
        for index, row in disease.iterrows():
            items = [row["DO qualifier"] + ": " + row["DO term"],
                    row["Based on orthology with (symbol)"],
                    row["Reference (FBrf ID)"]]
            items = [str(i) for i in items]
            r = html.Tr([html.Td(items[0]),
                         html.Td(items[1]),
                         html.Td(html.A(items[2], 
                                        href="http://flybase.org/reports/" 
                                        + items[2] + ".html",
                                        target="_blank"))
                ])
            ret.append(r)
    return ret


#### Processing page callbacks ####

@app.callback(
    Output('upload_raw_data_success_output', 'children'),
    [Input('upload_raw_data', 'contents')],
    [State('upload_raw_data', 'filename'),
     State('session-id', 'children')])
def parse_uploaded_10X(contents, filename, session_ID):
    if (filename is None):
        return dash.no_update

    print("[STATUS] parsing raw 10X data upload")
    if not ("zip" in filename):
        return "Uploaded file must be a flat .zip file containing the 10X directory contents"

    save_dir = helper_functions.save_analysis_path + "/" + session_ID + "/"
    if not (os.path.isdir(save_dir)):
        os.mkdir(save_dir)

    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)

    data = zf.ZipFile(io.BytesIO(decoded), mode="r")
    data.extractall(path=save_dir)

    return "Uploaded data parsed successfully"

@app.callback(
    Output('min_max_genes_slider_output_container', 'children'),
    [Input('min_max_genes_slider', 'value')]
)
def update_min_max_genes_output(value):
    print("[STATUS] [min_genes, max_genes] updated to " + str(value))
    return ("[min_genes, max_genes] = " + str(value))

@app.callback(
    Output('min_cells_slider_output_container', 'children'),
    [Input('min_cells_slider', 'value')]
)
def update_min_cells_output(value):
    print("[STATUS] min_cells updated to " + str(value))
    return ("min_cells = " + str(value))

@app.callback(
    Output('n_top_genes_slider_output_container', 'children'),
    [Input('n_top_genes_slider', 'value')]
)
def update_n_top_genes_output(value):
    print("[STATUS] n_top_genes updated to " + str(value))
    return ("# highly variable genes = " + str(value))

@app.callback(
    Output('n_neighbors_slider_output_container', 'children'),
    [Input('n_neighbors_slider', 'value')]
)
def update_n_neighbors_output(value):
    print("[STATUS] n_neighbors updated to " + str(value))
    return ("n_neighbors = " + str(value))

@app.callback(
    Output('clustering_resolution_slider_output_container', 'children'),
    [Input('clustering_resolution_slider', 'value')]
)
def update_clustering_resolution_output(value):
    print("[STATUS] clustering resolution updated to " + str(value))
    return "clustering resolution = " + str(value)

@app.callback(
    Output("processing_UMAP_plot", "figure"),
    [Input("refresh_all_button", "n_clicks"),
     Input("refresh_projection_button", "n_clicks"),
     Input("refresh_clustering_button", "n_clicks"),
     Input("refresh_pseudotime_button", "n_clicks"),
     Input("processing_UMAP_dropdown", "value")],
    [State("session-id", "children"),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value"),
     State("min_max_genes_slider", "value"),
     State("min_cells_slider", "value"),
     State("n_top_genes_slider", "value"),
     State("processing_UMAP_plot", "selectedData")]
)
def refresh_processing_UMAP(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, pt_btn_clicks, processing_plot_type, 
                            session_ID, n_neighbors, resolution, min_max_genes,
                            min_cells, n_top_genes, selected_cells,
                            adata=None, data_dir=None, target_sum=1e6, 
                            flavor="cell_ranger", n_comps=50, random_state=0):

    print("[STATUS] refreshing processing UMAP plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]


    if  (button_id == "refresh_clustering_button"):
        if (clust_btn_clicks in [None, 0]):
            return dash.no_update
        adata = cache_adata(session_ID)
        adata = do_clustering(session_ID, adata, resolution=resolution)
    
    elif(button_id == "refresh_projection_button"):
        if (proj_btn_clicks in [None, 0]):
            return dash.no_update
        adata = cache_adata(session_ID)   
        adata = do_neighborhood_graph(session_ID, adata, 
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state)
        adata = do_UMAP(session_ID, adata, random_state=random_state)
    
    elif(button_id == "refresh_all_button"):
        print("[DEBUG] refresh_all_button clicked")
        if (all_btn_clicks in [None, 0]):
            return dash.no_update
        print("[STATUS] refreshing everything")   
        adata = preprocess_data(session_ID, data_dir, min_cells=min_cells,
                                min_genes=min_max_genes[0], 
                                max_genes=min_max_genes[1], 
                                target_sum=target_sum, flavor=flavor, 
                                n_top_genes=n_top_genes)
        adata, = do_PCA(session_ID, adata, n_comps=n_comps, 
                        random_state=random_state),
        adata = do_neighborhood_graph(session_ID, adata, 
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state)
        adata = do_UMAP(session_ID, adata, random_state=random_state)
        adata = do_clustering(session_ID, adata, resolution=resolution)

    elif(button_id == "refresh_pseudotime_button"):
        print("[DEBUG] refresh_pseudotime_button clicked")
        if (pt_btn_clicks in [None, 0]):
            return dash.no_update
        
        # get the starter cell from the processing plot and redo the pseudotime based on that
        # if there are multiple selected, take the last one
        if not (selected_cells is None):
            adata = cache_adata(session_ID)
            pt_selected_cell_IDs = get_cell_intersection(session_ID, adata, [selected_cells])

        
        if ((pt_selected_cell_IDs is None) or (len(pt_selected_cell_IDs) != 1)):
            print("[ERROR] please select exactly 1 cell to use as a pseudotime starter cell"
                + "\nUsing first cell in set")
        starter_cell_ID = pt_selected_cell_IDs.pop()
        
        adata = do_pseudotime(session_ID, adata, starter_cell_ID)


    # if it's a dropdown menu update - load adata
    elif(button_id == "processing_UMAP_dropdown"):
        if (processing_plot_type in [0, "", None, []]):
            return dash.no_update
        else:
            adata = cache_adata(session_ID)

        # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return dash.no_update
    
    # update the plot
    if (processing_plot_type in [0, "", None, []]):
        return dash.no_update
    print("[STATUS] updating plot by: " + str(processing_plot_type))
    adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))

    if not (selected_cells is None):
        cells_to_highlight = get_cell_intersection(session_ID, adata, [selected_cells])
    else:
        cells_to_highlight = []
    if (processing_plot_type == "leiden_n"):
        return plotting.plot_UMAP(adata, "leiden_n", cells_to_highlight)
    elif (processing_plot_type == "pseudotime"):
        return plotting.plot_pseudotime_UMAP(adata, "pseudotime")
    elif (processing_plot_type == "differentiation potential"):
        return plotting.plot_pseudotime_UMAP(adata, "differentiation_potential")
    elif (processing_plot_type == "# UMIs (log1p)"):
        return plotting.plot_expression_UMAP(adata, "log1p_total_counts")
    elif (processing_plot_type == "# unique genes"):
        return plotting.plot_expression_UMAP(adata, "n_genes")

@app.callback(
    Output("processing_QC_plot", "figure"),
    [Input("processing_QC_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_violin_QC_plot(selected_QC, session_ID):
    if (selected_QC in [None, 0, []]):
        return dash.no_update

    print("[STATUS] updating violin gene plot")
    adata = cache_adata(session_ID)

    # plot function expects list of factors/genes, but for QC
    # we will only show one at a time - pass a list anyways though
    return plotting.plot_expression_violin(adata, [selected_QC])