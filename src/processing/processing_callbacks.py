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
    [Output(f"{x}-collapse", "is_open") for x in ["upload", "QC", "projection", "clustering"]],
    [Input(f"{x}-collapse-button", "n_clicks") for x in ["upload", "QC", "projection", "clustering"]],
    [State(f"{x}-collapse", "is_open") for x in ["upload", "QC", "projection", "clustering"]],
)
def toggle_procssing_accordion(n1, n2, n3, n4, is_open1, is_open2, is_open3, is_open4):
    print("[DEBUG] accordion triggered")
    ctx = dash.callback_context

    if not ctx.triggered:
        return dash.no_update
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    print("[DEBUG] button triggered: " + str(button_id))
    if button_id == "upload-collapse-button" and n1:
        return not is_open1, False, False, False
    elif button_id == "QC-collapse-button" and n2:
        return False, not is_open2, False, False
    elif button_id == "projection-collapse-button" and n3:
        return False, False, not is_open3, False
    elif button_id == "clustering-collapse-button" and n4:
        return False, False, False, not is_open4
    return False, False, False, False

@app.callback(
    Output('upload_raw_data_success_output', 'children'),
    [Input('upload_raw_data', 'contents')],
    [State('upload_raw_data', 'filename'),
     State('session-id', 'children')])
def parse_uploaded_data(contents, filename, session_ID):
    default_return = dash.no_update
    if (filename is None):
        return default_return

    n_steps = 3
    print("[STATUS] parsing data upload")
    
    content_type, content_string = contents.split(',')

    if (".h5ad" in filename):
        decoded = base64.b64decode(content_string)
        save_dir = save_analysis_path + "/" + str(session_ID) + "/"
        if not (os.path.isdir(save_dir)):
            os.makedirs(save_dir)
        with open(save_dir + "adata_cache.h5ad", "wb") as f:
            f.write(decoded)
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        adata = sc.read_h5ad(save_dir + "adata_cache.h5ad")
        #adata = cache_adata(session_ID, adata)
        adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
        adata.var_names_make_unique()

        state = {"filename": str(filename),
                 "# cells/obs": len(adata.obs.index),
                 "# genes/var": len(adata.var.index)}
        cache_state(session_ID, state)

        cache_adata(session_ID, adata)
        cache_progress(session_ID, progress=int(2/n_steps * 100))

        gene_list = adata.var.index.tolist()
        gene_list = [str(x) for x in gene_list]
        gene_list = list(sorted(gene_list, key=str.lower))
        cache_gene_list(session_ID, gene_list)
        cache_progress(session_ID, progress=int(3/n_steps * 100))
        cache_history(session_ID, history="Anndata object: " + str(filename) + " loaded successfully")

        return "Anndata object uploaded successfully"

    if (".zip" in filename):
        save_dir = save_analysis_path + "/" + str(session_ID) + "/raw_data/"
        if not (os.path.isdir(save_dir)):
            os.makedirs(save_dir)

        decoded = base64.b64decode(content_string)
        data = zf.ZipFile(io.BytesIO(decoded), mode="r")
        data.extractall(path=save_dir)
        cache_progress(session_ID, progress=int(1/n_steps * 100))


        adata = generate_adata_from_10X(session_ID)
        if (adata is None):
            return default_return
        adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
        adata.var_names_make_unique()
        
        state = {"filename": str(filename),
                 "# cells/obs": len(adata.obs.index),
                 "# genes/var": len(adata.var.index)}
        cache_state(session_ID, state)

        cache_adata(session_ID, adata)
        cache_progress(session_ID, progress=int(2/n_steps * 100))


        gene_list = adata.var.index.tolist()
        gene_list = [str(x) for x in gene_list]
        gene_list = list(sorted(gene_list, key=str.lower))
        cache_gene_list(session_ID, gene_list)
        cache_progress(session_ID, progress=int(3/n_steps * 100))
        cache_history(session_ID, history="10X data: " + str(filename) + " loaded successfully")

        return "Raw 10X data uploaded successfully"

    return "Uploaded file not recognized. Upload an anndata object in h5ad format or zipped 10X ouput data."

@app.callback(
    Output('load_selected_dataset_success_output', 'children'),
    [Input('processing_load_dataset_button', 'n_clicks')],
    [State('processing_dataset_dropdown', 'value'),
     State('session-id', 'children')])
def parse_selected_dataset(btn_clicks, dataset_key, session_ID):
    default_return = ""
    print("[STATUS] loading selected dataset")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if  (button_id == "processing_load_dataset_button"):
        if (btn_clicks in [None, 0]):
            return default_return

    if (dataset_key in [0, "", None, []]):
        return default_return

    n_steps = 2
    adata = load_selected_dataset(session_ID, dataset_key)
    if (adata is None):
        print("[ERROR] selected dataset not found; key: " + str(dataset_key))
        return default_return
    else:
        adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
        cache_adata(session_ID, adata)
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        gene_list = adata.var.index.tolist()
        gene_list = [str(x) for x in gene_list]
        gene_list = list(sorted(gene_list, key=str.lower))
        cache_gene_list(session_ID, gene_list)
        cache_progress(session_ID, progress=int(2/n_steps * 100))
        cache_history(session_ID, history="Curated dataset loaded")

        return "Dataset loaded"

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
    [Output("processing_UMAP_plot", "figure"),
     Output("refresh_all_status", "children")],
    [Input("refresh_all_button", "n_clicks"),
     Input("refresh_projection_button", "n_clicks"),
     Input("refresh_clustering_button", "n_clicks"),
     Input("processing_UMAP_dropdown", "value"),
     Input("neighbors_method_radio", "value"),
     Input("n_dims_processing_radio", "value")],
    [State("session-id", "children"),
     State("n_neighbors_slider", "value"),
     State("clustering_resolution_slider", "value"),
     State("min_max_genes_slider", "value"),
     State("min_cells_slider", "value"),
     State("n_top_genes_slider", "value")]
)
def refresh_processing_UMAP(all_btn_clicks, proj_btn_clicks,
                            clust_btn_clicks, processing_plot_type,
                            neighborhood_method, n_dim_proj_plot, 
                            session_ID, n_neighbors, resolution, 
                            min_max_genes, min_cells, n_top_genes,
                            adata=None, target_sum=1e6, 
                            flavor="cell_ranger", n_comps=50, random_state=0):
    
    default_return = [dash.no_update, dash.no_update]
    print("[STATUS] refreshing processing UMAP plot")
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
        print("[DEBUG] refresh_all_button clicked")
        if (all_btn_clicks in [None, 0]):
            return default_return

        n_steps = 6
        adata = cache_adata(session_ID)
        if (adata is None):
            return default_return
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        print("[STATUS] refreshing everything")   
        adata = preprocess_data(session_ID, adata, min_cells=min_cells,
                                min_genes=min_max_genes[0], 
                                max_genes=min_max_genes[1], 
                                target_sum=target_sum, flavor=flavor, 
                                n_top_genes=n_top_genes)
        cache_progress(session_ID, progress=int(2/n_steps * 100))

        adata, = do_PCA(session_ID, adata, n_comps=n_comps, 
                        random_state=random_state),
        cache_progress(session_ID, progress=int(3/n_steps * 100))

        adata = do_neighborhood_graph(session_ID, adata, neighborhood_method,
                                      n_neighbors=n_neighbors, 
                                      random_state=random_state)
        cache_progress(session_ID, progress=int(4/n_steps * 100))

        adata = do_UMAP(session_ID, adata, random_state=random_state)
        cache_progress(session_ID, progress=int(5/n_steps * 100))

        adata = do_clustering(session_ID, adata, resolution=resolution)
        cache_progress(session_ID, progress=int(6/n_steps * 100))
        
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

    print("[STATUS] updating plot by: " + str(processing_plot_type))
    #adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))

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

    print("[STATUS] updating violin gene plot")
    if (adata_cache_exists(session_ID) is False):
        return default_return

    # plot function expects list of factors/genes, but for QC
    # we will only show one at a time here - list is required though
    return plot_expression_violin(session_ID, [selected_QC], show_points = False)