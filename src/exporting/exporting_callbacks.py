import os
import io

from flask import send_file

from app import app
from helper_functions import *


#### Export analysis page callbacks ####
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
                             mimetype="application/octet-stream"
                    )
    else:
        print("[ERROR] file " + f + " not available for download")
        return ("Error - file not generated yet. Go back in your browser," 
              + "upload your raw data, and perform QC before downloading.")