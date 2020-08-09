import dash
from dash.dependencies import Input, Output, State

import shortuuid

import os
import io

from flask import send_file
from flask_security import current_user

from app import app
from helper_functions import *

'''
@app.callback(
    [Output("exporting_button", "disabled"),
     Output("exporting_link", "href")],
    [Input("main_tabs", "active_tab")],
    [State("session-id", "children")]
)
def enabled_export_button(active_tab, session_ID):
    default_return = [True, dash.no_update] #disabled
    
    if (active_tab != "exporting_tab"):
        return default_return
    
    if not (adata_cache_exists(session_ID)):
        return default_return
    else:
        return [False, "/download/" + session_ID + "/h5ad"]
'''

@app.callback(
    [Output("exporting_obs_column_dropdown", "disabled"),
     Output("exporting_obs_value_dropdown", "disabled")],
    [Input("exporting_subset_radio", "value")]
)
def enabled_export_button(subset_cell_option):    
    if (subset_cell_option == "all"):
        return [True, True] #enabled
    else:
        return [False, False] #disabled

@app.callback(
    Output("exporting_obs_column_dropdown", "options"),
    [Input("main_tabs", "active_tab")],
    [State('session-id', 'children')]
)
def update_exporting_obs_column_dropdown(active_tab, session_ID):
    default_return = dash.no_update

    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (active_tab != "exporting_tab"):
        return default_return

    adata = cache_adata(session_ID)
    if (adata is None):
        return default_return

    a = adata.obs.select_dtypes(include=["category"])
    options = [
        {"label": str(x), "value": str(x)} for x in a.columns.to_list()
    ]
    return options


@app.callback(
    Output("exporting_obs_value_dropdown", "options"),
    [Input("exporting_obs_column_dropdown", "value")],
    [State("session-id", "children")]
)
def update_exporting_obs_value_dropdown(obs_column, session_ID):
    default_return = dash.no_update

    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (obs_column in [None, [], 0, ""]):
        return default_return

    adata = cache_adata(session_ID)
    if (adata is None):
        return default_return

    vals = adata.obs[obs_column].unique()
    options = [
        {"label": str(x), "value": str(x)} for x in list(vals)
    ]
    return options


@app.callback(
    [Output("exporting_obs_link_button", "disabled"),
     Output("exporting_obs_subset_link", "href")],
    [Input("exporting_obs_button", "n_clicks")],
    [State("session-id", "children"),
     State("exporting_obs_column_dropdown", "value"),
     State("exporting_obs_value_dropdown", "value"),
     State("exporting_subset_radio", "value")]
)
def obs_subset_anndata_h5ad(n_clicks, session_ID,
                            obs_column, obs_values, subset_cell_option):
    default_return = [True, dash.no_update]

    if (n_clicks in [0, "", None, []]):
        return default_return

    if not (subset_cell_option == "all"):
        if (obs_column in [0, "", None, []]):
            return default_return

        if (obs_values in [0, "", None, []]):
            return default_return

    adata = cache_adata(session_ID)
    if (adata is None):
        return default_return

    if (subset_cell_option == "all"):
        out_adata = adata
    else:
        out_adata = adata[adata.obs[obs_column].isin(obs_values)]
    
    f = save_analysis_path + str(session_ID) + "/adata_export.h5ad"
    if (os.path.isfile(f)):
        os.remove(f)

    out_adata.write(f)
    ret = [False, "/download/" + session_ID + "/export_h5ad"]

    return ret

@app.callback(
    [Output("save_dataset_name", "valid"),
     Output("save_dataset_name", "invalid"),
     Output("save_dataset_output", "children")],
    [Input("save_dataset_button", "n_clicks")],
    [State("session-id", "children"),
     State("save_dataset_name", "value"),
     State("exporting_obs_column_dropdown", "value"),
     State("exporting_obs_value_dropdown", "value"),
     State("exporting_subset_radio", "value")]
)
def save_anndata_dataset(n_clicks, session_ID, save_dataset_name,
                         obs_column, obs_values, subset_cell_option):
    default_return = [dash.no_update, dash.no_update, dash.no_update]
    failed_return = [False, True, dash.no_update]
    successful_return = [True, False, dash.no_update]

    if (n_clicks in [0, "", None, []]):
        return default_return

    if not (subset_cell_option == "all"):
        if (obs_column in [0, "", None, []]):
            return failed_return

        if (obs_values in [0, "", None, []]):
            return failed_return

    adata = cache_adata(session_ID)
    if (adata is None):
        return failed_return

    if (subset_cell_option == "all"):
        out_adata = adata
    else:
        out_adata = adata[adata.obs[obs_column].isin(obs_values)]
    
    store_dir = user_dataset_path + str(current_user.email) + "/"
    store_name = save_dataset_name
    print("[DEBUG] attempting to save at: " 
         + str(store_dir) + str(store_name))

    cache_adata(session_ID, adata=out_adata, 
                store_dir=store_dir, store_name=store_name)
    ret = successful_return
    return ret

@app.server.route('/download/<path:path>/export_h5ad')
def serve_subset_anndata_h5ad(path):
    f = save_analysis_path + path + "/adata_export.h5ad"
    if (os.path.isfile(f)):
        with open(f, "rb") as b:
            return send_file(io.BytesIO(b.read()), 
                      as_attachment=True,
                      attachment_filename="adata_export.h5ad",
                      mimetype="application/octet-stream", 
                      cache_timeout=3
                    )
    else:
        print("[ERROR] issue with subsetting; file not found after saving subset.")
        return "[ERROR] issue with subsetting; file not found after saving subset."

'''
@app.server.route('/download/<path:path>/h5ad')
def serve_anndata_h5ad(path):
    print("[STATUS] preparing to serve anndata in h5ad format")
    f = save_analysis_path + path + "/adata_cache.h5ad"
    if (os.path.isfile(f)):
        print("[DEBUG] file " + f + " found - serving")
        with open(f, "rb") as b:
            return send_file(io.BytesIO(b.read()), 
                             as_attachment=True,
                             attachment_filename="adata.h5ad",
                             mimetype="application/octet-stream", cache_timeout=3
                    )
    else:
        print("[ERROR] file " + f + " not available for download")
        return ("Error - file not generated yet. Go back in your browser," 
             + "upload your raw data, and perform QC before downloading.")
'''