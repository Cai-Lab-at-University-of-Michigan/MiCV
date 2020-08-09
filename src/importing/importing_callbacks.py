import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from flask_security import current_user

import os
import shutil
import zipfile as zf
import base64
import io

from helper_functions import *
from status.status_functions import *
from app import app

from . importing_functions import *

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
    
    content_type, content_string = contents.split(',')

    if (".h5ad" in filename):
        decoded = base64.b64decode(content_string)
        save_dir = save_analysis_path + "/" + str(session_ID) + "/"
        if not (os.path.isdir(save_dir)):
            os.makedirs(save_dir)
        print("[DEBUG] attempting to save uploaded data")
        with open(save_dir + "adata_cache.h5ad", "wb") as f:
            f.write(decoded)
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        adata = sc.read_h5ad(save_dir + "adata_cache.h5ad")
        #adata = cache_adata(session_ID, adata)
        adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
        adata.var_names_make_unique()

        # figure out how many total UMIs are here
        if ("total_counts" in adata.obs):
            n_counts = int(np.sum(adata.obs["total_counts"]))
        else: # assume uploaded values in X represent UMI counts
            n_counts = int(np.sum(adata.X))

        state = {"filename": str(filename),
                 "# cells/obs": len(adata.obs.index),
                 "# genes/var": len(adata.var.index),
                 "# counts": n_counts}
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
                 "# genes/var": len(adata.var.index),
                 "# counts": np.sum(adata.X)} # assume uploaded values in X are UMI counts
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
    [Input('importing_load_dataset_button', 'n_clicks')],
    [State('importing_dataset_dropdown', 'value'),
     State('session-id', 'children')])
def parse_selected_dataset(btn_clicks, dataset_key, session_ID):
    default_return = ""
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if  (button_id == "importing_load_dataset_button"):
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
        cache_adata(session_ID, adata.obs, group="obs")
        cache_progress(session_ID, progress=int(1/n_steps * 100))

        gene_list = adata.var.index.tolist()
        gene_list = [str(x) for x in gene_list]
        gene_list = list(sorted(gene_list, key=str.lower))
        cache_gene_list(session_ID, gene_list)
        cache_progress(session_ID, progress=int(2/n_steps * 100))
        cache_history(session_ID, history="Curated dataset loaded")

        return "Dataset loaded"

@app.callback(
    Output('importing_user_dataset_success_output', 'children'),
    [Input('importing_user_dataset_load_button', 'n_clicks'),
     Input('importing_user_dataset_delete_button', 'n_clicks')],
    [State('importing_user_dataset_list_radio', 'value'),
     State('session-id', 'children')])
def access_user_dataset(load_btn_clicks, delete_btn_clicks,
                        dataset_name, session_ID):
    default_return = dash.no_update
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if (dataset_name in [0, "", None, []]):
        return default_return
    else:
        dataset_dir = (user_dataset_path 
               + str(current_user.email) + "/")

    if (button_id == "importing_user_dataset_load_button"):
        if (load_btn_clicks in [None, 0]):
            return default_return
        else:
            # load the data
            n_steps = 2
            adata = cache_adata(session_ID, store_dir=dataset_dir, store_name=dataset_name)
            if (adata is None):
                print("[ERROR] user dataset not found; name: " + str(dataset_name))
                return default_return
            else:
                cache_adata(session_ID, adata)
                cache_progress(session_ID, progress=int(1/n_steps * 100))

                gene_list = adata.var.index.tolist()
                gene_list = [str(x) for x in gene_list]
                gene_list = list(sorted(gene_list, key=str.lower))
                cache_gene_list(session_ID, gene_list)
                cache_state(session_ID, key="filename", val=str(dataset_name))
                cache_progress(session_ID, progress=int(2/n_steps * 100))
                cache_history(session_ID, history="User dataset loaded")

                return "Dataset loaded"
    if (button_id == "importing_user_dataset_delete_button"):
        if (delete_btn_clicks in [None, 0]):
            return default_return
        else:
            if (os.path.isdir(dataset_dir) is True):
                shutil.rmtree(dataset_dir, ignore_errors=True)
            return "Dataset deleted"

@app.callback(
    [Output("importing_user_dataset_list_radio", "options"),
     Output("importing_user_dataset_load_button", "disabled"),
     Output("importing_user_dataset_delete_button", "disabled")],
    [Input("status_interval", "n_intervals")])
def update_user_dataset_list(n):
    default_return = [[{"label": "You have no datasets saved", 
                        "value": "no_datasets_saved"}],
                       True, True]
    if (current_user is None):
        raise PreventUpdate
    try:
        # lookup the user's data store by email
        dataset_list = get_dataset_list(current_user.email)
    except:
        return default_return

    if (dataset_list in [[], None, 0, ""]):
        return default_return
    else:
        options=[{"label": d[0], "value": d[0]} for d in dataset_list]
        return [options, False, False]