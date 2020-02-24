import dash
from dash.dependencies import Input, Output, State

import shortuuid

import os
import io

from flask import send_file

from app import app
from helper_functions import *

@app.callback(
     Output("exporting_obs_value_dropdown", "options"),
    [Input("exporting_obs_column_dropdown", "value")],
    [State("session-id", "children")]
)
def update_exporting_obs_value_dropdown(obs_column, session_ID):
    default_return = dash.no_update

    print("[STATUS] refreshing marker gene UMAP dropdown options")
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
    Output("exporting_obs_link_button", "disabled"),
    [Input("exporting_obs_button", "n_clicks")],
    [State("session-id", "children"),
     State("exporting_obs_column_dropdown", "value"),
     State("exporting_obs_value_dropdown", "value")]
)
def obs_subset_anndatata_h5ad(n_clicks, session_ID,
                              obs_column, obs_values):
    default_return = True

    if (n_clicks in [0, "", None, []]):
        return default_return

    if (obs_column in [0, "", None, []]):
        return default_return

    if (obs_values in [0, "", None, []]):
        return default_return

    adata = cache_adata(session_ID)
    out_adata = adata[adata.obs[obs_column].isin(obs_values)]
    f = save_analysis_path + "adata_subset.h5ad"
    if (os.path.isfile(f)):
        os.remove(f)

    out_adata.write(f)

    return False
    
@app.server.route('/MiCV/download/subset_h5ad')
def serve_subset_anndata_h5ad():
    f = save_analysis_path + "adata_subset.h5ad"
    if (os.path.isfile(f)):
        with open(f, "rb") as b:
            return send_file(io.BytesIO(b.read()), 
                      as_attachment=True,
                      attachment_filename="adata_subset.h5ad",
                      mimetype="application/octet-stream", cache_timeout=3
                    )
    else:
        print("[ERROR] issue with subsetting; file not found after saving subset.")
        return "[ERROR] issue with subsetting; file not found after saving subset."


@app.server.route('/MiCV/download/h5ad')
def serve_anndata_h5ad():
    print("[STATUS] preparing to serve anndata in h5ad format")
    f = save_analysis_path + "adata_cache.h5ad"
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