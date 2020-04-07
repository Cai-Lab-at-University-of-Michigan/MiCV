import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html

from pseudotime.pseudotime_functions import calculate_gene_trends

from helper_functions import *
from plotting.plotting_functions import *
from app import app

#### Annotation page callbacks ####

@app.callback(
    [Output("clustering_UMAP_plot", "figure"),
     Output("total_cell_count", "value")],
    [Input("main_tabs", "active_tab"),
     Input("define_cluster_button", "n_clicks"),
     Input("clustering_dropdown", "value"),
     Input("pseudotime_UMAP_plot", "selectedData"),
     Input("expression_UMAP_plot", "selectedData"),
     Input("violin_gene_plot", "selectedData"),
     Input("pseudotime_gene_plot", "selectedData")],
    [State("session-id", "children"),
     State("clustering_UMAP_plot", "selectedData")]
)
def refresh_clustering_plot(active_tab, 
                            cluster_btn_clicks, clustering_plot_type, pt_selected,
                            expr_selected, violin_selected, pt_gene_selected, 
                            session_ID, clust_selected,  
                            adata=None, data_dir=None):

    default_return = [dash.no_update, dash.no_update]
    print("[STATUS] refreshing clustering plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (clustering_plot_type in ["", 0, None, []]):
        return default_return

    if(button_id == "annotation_tab"):
        if (active_tab != "annotation_tab"):
            return default_return
        else:
            obs = cache_adata(session_ID, group="obs")
            var = cache_adata(session_ID, group="var")
            if ((obs is None) or (var is None)):
                return default_return
            else:
                obs["leiden_n"] = pd.to_numeric(obs["leiden"])
                obs["cell_ID"] = obs.index
                obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(obs.index))))
                cache_adata(session_ID, obs, group="obs")

                gene_trends = cache_gene_trends(session_ID)  

                gene_list = var.index.tolist()
                gene_list = [str(x) for x in gene_list]
                gene_list = list(sorted(gene_list, key=str.lower))
                cache_gene_list(session_ID, gene_list)

            if (clustering_plot_type in [0, "", None, []]):
                return default_return
    
    elif(button_id == "define_cluster_button"):
        if (cluster_btn_clicks in [0, None, "", []]):
            return default_return

        obs = cache_adata(session_ID, group="obs")
        clustering_group = clustering_plot_type
        if not (clustering_group == "leiden"):
            if not (clustering_group in obs):
                obs[clustering_group] = "0"
            obs[clustering_group] = obs[clustering_group].astype(str)
            # create a new unique cluster ID for this cluster in this clustering_group
            previous_cluster_IDs = obs[clustering_group].unique()
            new_cluster_ID = str(len(previous_cluster_IDs))
            i = 0
            while (new_cluster_ID in previous_cluster_IDs):
                new_cluster_ID = str(len(previous_cluster_IDs) + i) 
                i += 1

            # for violin plot selection, take the intersection of all points selected
            # i.e. if the user selected cells from multiple expression violin plots,
            # take the intersection of those cells (logical AND)
            violin_selected = get_violin_intersection(session_ID, violin_selected) 
            
            pt_min, pt_max = get_pseudotime_min_max(session_ID, pt_gene_selected)

            # add the selected cells to this cluster
            selected_cells = get_cell_intersection(session_ID, 
                             [clust_selected, pt_selected,
                             expr_selected, violin_selected],
                             pt_min, pt_max)
    
            obs.loc[selected_cells, clustering_group] = new_cluster_ID
            obs[clustering_group] = obs[clustering_group].astype('category')
            cache_adata(session_ID, obs, group="obs")

    elif(button_id in ["Pseudotime_UMAP_plot",
                       "Expression_UMAP_plot",
                       "Violin_gene_plot",
                       "Pseudotime_gene_plot"]):
        if (pt_selected is None
        and expr_selected is None
        and violin_selected is None
        and pt_gene_selected is None):
            default_return


    # if it's a dropdown menu update - load adata
    elif(button_id == "clustering_dropdown"):
        if (clustering_plot_type in [0, "", None, []]):
            default_return

    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        default_return


    obs = cache_adata(session_ID, group="obs")
    if ((obs is None)
    or not (clustering_plot_type in obs)):
        return default_return

    # figure out which cells need to be selected, based on other graphs
    violin_selected = get_violin_intersection(session_ID, violin_selected)
    pt_min, pt_max = get_pseudotime_min_max(session_ID, pt_gene_selected)
    selected_cell_intersection = get_cell_intersection(session_ID,
                                                        [clust_selected, 
                                                         pt_selected,
                                                         expr_selected,
                                                         violin_selected],
                                                         pt_min, pt_max)
    selected_points = []
    for c in selected_cell_intersection:
        selected_points.append((obs).index.get_loc(c))
    if (len(selected_points) == 0):
        selected_points = range(0, len((obs.loc[:, "cell_ID"]).index.to_list()))
    
    # update the plot
    return plot_UMAP(session_ID, clustering_plot_type, selected_cell_intersection), len(obs.index.tolist())

@app.callback(
    Output("clustering_UMAP_count", "children"),
    [Input("clustering_UMAP_plot", "selectedData")],
    [State("session-id", "children"),
     State("total_cell_count", "value")]
)
def refresh_UMAP_clustering_count(selected_cells, session_ID, n_total_cells):
    if (selected_cells in ["", 0, [], None]):
        return dash.no_update

    if (n_total_cells in [0, None, "", []]):
        n_total_cells = 1
    n_cells_selected = len(selected_cells["points"])
    return ("# cells selected: " + str(n_cells_selected) + " | % total: "
             + str(round(100.0 * n_cells_selected/n_total_cells, 2)))


@app.callback(
    Output("pseudotime_UMAP_plot", "figure"),
    [Input("pseudotime_dropdown", "value")],
    [State('session-id', 'children'),
     State("main_tabs", "active_tab"),
     State("n_dims_proj_expression_radio", "value")]
)
def refresh_pseudotime_plot(pt_plot_type, session_ID, active_tab, 
                            n_dims_proj, adata=None, data_dir=None):

    print("[STATUS] refreshing pseudotime plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (button_id == "pseudotime_dropdown"):
        if (active_tab != "annotation_tab"):
            return dash.no_update

    # do nothing if no buttons pressed
    elif((button_id == "not_triggered")
      or (pt_plot_type in [0, "", None, []])):
        return dash.no_update
    
    # get the adata object from the cache
    obs = cache_adata(session_ID, group="obs")

    if ((obs is None)
    or not (pt_plot_type in obs)):
        return dash.no_update
    
    # regardless of what updates were requested - update the plot
    print("[STATUS] updating pseudotime plot")
    return plot_pseudotime_UMAP(session_ID, pt_plot_type)

@app.callback(
    Output("pseudotime_UMAP_count", "children"),
    [Input("pseudotime_UMAP_plot", "selectedData")],
    [State("session-id", "children"),
     State("total_cell_count", "value")]
)
def refresh_UMAP_pseudotime_count(selected_cells, session_ID, n_total_cells):
    if (selected_cells in ["", 0, [], None]):
        return dash.no_update
    
    if (n_total_cells in [0, None, "", []]):
        n_total_cells = 1

    n_cells_selected = len(selected_cells["points"])
    return ("# cells selected: " + str(n_cells_selected) + " | % total: "
             + str(round(100.0 * n_cells_selected/n_total_cells, 2)))


@app.callback(
    Output("single_gene_dropdown", "options"),
    [Input("main_tabs", "active_tab")],
    [State("session-id", "children")]
)
def update_single_gene_dropdown(active_tab, session_ID):
    print("[DEBUG] active_tab: " + str(active_tab))
    if (active_tab != "annotation_tab"):
        return dash.no_update
    
    gene_list = cache_gene_list(session_ID)
    if (gene_list in [None, 0, "", []]):
        return dash.no_update

    options = [{"label": i, "value": i} for i in gene_list]
    return options

@app.callback(
    Output("mixed_gene_dropdown", "options"),
    [Input("main_tabs", "active_tab")],
    [State("session-id", "children")]
)
def update_mixed_gene_dropdown(active_tab, session_ID):
    if (active_tab != "annotation_tab"):
        return dash.no_update
    gene_list = cache_gene_list(session_ID)
    if (gene_list in [None, 0, "", []]):
        return dash.no_update
    
    options = [{"label": i, "value": i} for i in gene_list]
    return options

@app.callback(
    [Output("expression_UMAP_plot", "figure"),
     Output("mixed_expression_legend_image", "hidden")],
    [Input("single_gene_dropdown", "value"),
     Input("mixed_gene_dropdown", "value"),
     Input("single_gene_expression_radio", "value"),
     Input("n_dims_proj_expression_radio", "value")],
    [State('session-id', 'children')]
)
def refresh_expression_UMAP_plot(selected_gene, selected_mixed_genes,
                                 multi, n_dims_proj, session_ID):
    default_return = [dash.no_update, dash.no_update]

    if (multi == "standard"):
        if (selected_gene in [None, 0, []]):
            return default_return
        else:
            plot_these_genes = selected_gene
            hidden = True

    else:
        if (selected_mixed_genes in [None, 0, []]):
            return default_return
        else:
            plot_these_genes = selected_mixed_genes
            hidden = False


    if (adata_cache_exists(session_ID) is False):
        return default_return

    print("[STATUS] updating expression UMAP plot")
    return [plot_expression_UMAP(session_ID, plot_these_genes, multi, n_dim=n_dims_proj), hidden]

@app.callback(
    Output("gene_UMAP_count", "children"),
    [Input("expression_UMAP_plot", "selectedData")],
    [State("session-id", "children"),
     State("total_cell_count", "value")]
)
def refresh_UMAP_gene_count(selected_cells, session_ID, n_total_cells):
    default_return = dash.no_update

    if (selected_cells in ["", 0, [], None]):
        return default_return

    n_cells_selected = len(selected_cells["points"])
    return ("# cells selected: " + str(n_cells_selected) + " | % total: "
             + str(round(100.0 * n_cells_selected/n_total_cells, 2)))


@app.callback(
    Output("multi_gene_dropdown", "options"),
    [Input("main_tabs", "active_tab")],
    [State("session-id", "children")]
)
def update_multi_gene_dropdown(active_tab, session_ID):
    if (active_tab != "annotation_tab"):
        return dash.no_update
    
    gene_list = cache_gene_list(session_ID)
    if (gene_list in [None, 0, "", []]):
        return dash.no_update

    options = [{"label": i, "value": i} for i in gene_list]
    return options

@app.callback(
    [Output("pseudotime_gene_plot", "figure"),
     Output("violin_gene_plot", "figure")],
    [Input("multi_gene_dropdown", "value"),
     Input("pseudotime_gene_relative_radio", "value"),
     Input("pseudotime_gene_branch_dropdown", "value")],
    [State('session-id', 'children')]
)
def refresh_pseudotime_gene_plot(selected_genes, relative, branch_n, session_ID):
    default_return = [dash.no_update, dash.no_update]
    if (selected_genes in [None, 0, [], ""]):
        return default_return
    
    ret = default_return

    print("[STATUS] plotting gene pseudotime expression for: " + str(selected_genes))
    print("[DEBUG] branch_n: " + str(branch_n))
    if not(branch_n in ["", None, []]):
        #gene_trends = cache_gene_trends(session_ID)
        print("[STATUS] calculating gene pseudotime expression for: " + str(selected_genes))
        gene_trends = calculate_gene_trends(session_ID, selected_genes, branch_n)
        if not (gene_trends in ["", 0, None, []]):
            branch = list(gene_trends.keys())[0] #branch names (cell_IDs of terminal cells)
            ret[0] = plot_expression_trend(gene_trends, selected_genes, 
                                        selected_branch=branch, 
                                        relative=relative)
    
    print("[STATUS] updating violin gene plot")
    ret[1] = plot_expression_violin(session_ID, selected_genes)
    return ret


@app.callback(
    Output("gene_violin_count", "children"),
    [Input("violin_gene_plot", "selectedData")],
    [State("session-id", "children"),
     State("total_cell_count", "value")]
)
def refresh_violin_gene_count(selected_cells, session_ID, n_total_cells):
    if (selected_cells in ["", 0, [], None]):
        return dash.no_update
    
    if (n_total_cells in [0, None, "", []]):
        n_total_cells = 1

    violin_selected = get_violin_intersection(session_ID, selected_cells)
    n_cells_selected = len(violin_selected["points"])
    return ("# cells selected: " + str(n_cells_selected) + " | % total: "
             + str(round(100.0 * n_cells_selected/n_total_cells, 2)))


@app.callback(
    Output('gene_data_table', 'children'),
    [Input('single_gene_dropdown', 'value')],
    [State('session-id', 'children'),
     State('single_gene_expression_radio', 'value')]
)
def update_gene_data_table(selected_gene, session_ID, multi):
    print("[STATUS] getting data for first of selected_gene: " 
          + str(selected_gene))

    if (selected_gene in [None, "", 0, []]):
        return dash.no_update

    if (multi == "multiple"):
        selected_gene = selected_gene[0]

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