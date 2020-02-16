import pandas as pd
import palantir

from helper_functions import *

def do_pseudotime(session_ID, adata, starter_cell_ID=None):
    a = adata.to_df()

    print("[STATUS] computing pseudotime ...")
    pca_projections, var_r = palantir.utils.run_pca(a)

    dm_res = palantir.utils.run_diffusion_maps(pca_projections, knn=20)
    ms_data = palantir.utils.determine_multiscale_space(dm_res)

    # We're just going to use the UMAP projection instead
    tsne = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs.index)
    tsne.columns = ["x", "y"]

    imp_df = palantir.utils.run_magic_imputation(a, dm_res, n_steps=1)
    
    #start_cell = str(adata.obs.index[0])
    if (starter_cell_ID is None):
        start_cell = "TTGTTTGCAATTTCCT"
        print("[ERROR] no starter cell provided; using " + start_cell +
              " as a default (assuming original type-II data)")
    else:
        start_cell = starter_cell_ID
    print("[DEBUG] start cell is: " + str(start_cell))
    pr_res = palantir.core.run_palantir(ms_data, start_cell, 
                                        terminal_states=None, knn=20, 
                                        num_waypoints=500, 
                                        use_early_cell_as_start=True, 
                                        scale_components=False)
    adata.obs["pseudotime"] = pr_res.pseudotime[tsne.index]
    adata.obs["differentiation_potential"] = pr_res.entropy[tsne.index]

    genes = adata.var.index.tolist()
    #genes = ["CycE", "dap", "Hey", "nSyb"]
    #print("[DEBUG] genes: " + str(genes))
    #print("[DEBUG] d.imp_df.loc[:, genes]: " + str(d.imp_df.loc[:, genes]))
    print("[STATUS] computing all gene trends (this will take a while)")
    gene_trends = palantir.presults.compute_gene_trends(pr_res, 
                                                        imp_df.loc[:, genes])

    cache_adata(session_ID, adata)
    cache_gene_trends(session_ID, gene_trends)
    return adata