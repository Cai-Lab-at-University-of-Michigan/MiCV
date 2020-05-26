import pandas as pd
import numpy as np
from pygam import LinearGAM, ExpectileGAM, s, f
import palantir

from helper_functions import *
from status.status_functions import *

def do_pseudotime(session_ID, adata, starter_cell_ID=None):
    n_steps = 6
    adata_highly_variable = adata[:, adata.var["highly_variable"]]
    a = adata_highly_variable.to_df()
    cache_progress(session_ID, progress=int(1/n_steps * 100))
    print("[STATUS] computing pseudotime ...")

    pca_projections, var_r = palantir.utils.run_pca(a)
    cache_progress(session_ID, progress=int(2/n_steps * 100))
    cache_history(session_ID, history=("Performed PCA for pseudotime"))

    dm_res = palantir.utils.run_diffusion_maps(pca_projections, knn=20)
    cache_progress(session_ID, progress=int(3/n_steps * 100))

    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    cache_progress(session_ID, progress=int(4/n_steps * 100))
    cache_history(session_ID, history=("Generated diffusion maps and"
                                     + " multiscale space for pseudotime"))

    #imp_df = palantir.utils.run_magic_imputation(a, dm_res, n_steps=2)
    adata.layers["imputed"] = sc.external.pp.magic(adata, t=2, 
                                                   solver="approximate",
                                                   n_pca = 19, 
                                                   n_jobs=1, copy=True).X
    cache_progress(session_ID, progress=int(5/n_steps * 100))
    cache_history(session_ID, history=("Imputed data for pseudotime"))
    
    start_cell = starter_cell_ID
    # TODO: need to catch ValueError here in try-catch block
    try:
        pr_res = palantir.core.run_palantir(ms_data, start_cell, 
                                            terminal_states=None,
                                            num_waypoints=1200, 
                                            scale_components=True, n_jobs=1)
    except ValueError:
        cache_progress(session_ID, progress=100)
        cache_history(session_ID, history=("[ERROR] pseudotime failed to converge. Try again with"
                                         + " a different starter cell."))
        return adata

    cache_progress(session_ID, progress=int(5/n_steps * 100))
    cache_history(session_ID, history=("Calculated pseudotime trajectories with "
                                     + "starter cell ID: " + str(starter_cell_ID)))

    adata.obs["pseudotime"] = pr_res.pseudotime[adata.obs.index]
    adata.obs["differentiation_potential"] = pr_res.entropy[adata.obs.index]
    
    # drop old branch probabilities if they exist, then add new ones
    adata.obs.drop(list(adata.obs.filter(regex='pseudotime_branch_')), 
                   axis=1, inplace=True)
    for i, branch in enumerate(pr_res.branch_probs.columns):
        adata.obs["pseudotime_branch_" + str(i)] = pr_res.branch_probs.loc[adata.obs.index, branch]
        #adata.obs["in_pseudotime_branch_" + str(i)] = adata.obs["pseudotime_branch_" + str(i)] >= 0.7

    adata.uns["branch_probs"] = pr_res.branch_probs.to_dict()
    adata.uns["pseudotime"] = pr_res.pseudotime.to_dict()
    cache_adata(session_ID, adata)

    cache_progress(session_ID, progress=int(6/n_steps * 100))
    return adata

def calculate_gene_trends(session_ID, list_of_genes, branch_ID):
    n_steps = 2 + len(list_of_genes)
    
    uns = cache_adata(session_ID, group="uns")
    obs = cache_adata(session_ID, group="obs")
    cache_progress(session_ID, progress=int(1/n_steps * 100))

    if ("branch_probs" in uns.keys()):
        branch_probs = pd.DataFrame.from_dict(uns["branch_probs"])
    else:
        branch_probs = None

    pseudotime = obs["pseudotime"]
    cache_progress(session_ID, progress=int(2/n_steps * 100))

    if ((branch_ID == -1) or (branch_probs is None)):
        branch = "all branches"
        cells_in_branch = obs.index
    else:
        branch = list(branch_probs.columns)[branch_ID]
        cells_in_branch = obs[obs["pseudotime_branch_" + str(branch_ID)] > 0.2].index
    print("[DEBUG] branch: " + str(branch))

    '''
    gene_trends = palantir.presults.compute_gene_trends(pr_res, 
                                                        imp_df.loc[:, genes],
                                                        lineages = [branch],
                                                        n_jobs=1)
    '''
    X_train = pseudotime.to_numpy()
    if ((branch != "all branches") and not (branch_probs is None)):
        weights = branch_probs[branch].to_numpy()
    else:
        weights = np.ones_like(X_train)

    X_train = np.reshape(X_train, (len(X_train), 1))
    weights = np.reshape(weights, (len(weights), 1))

    X_plot = np.linspace(np.min(obs["pseudotime"][cells_in_branch]), 
                         np.max(obs["pseudotime"][cells_in_branch]), 
                         250)
    

    gene_trends = pd.DataFrame()
    gene_trends["pseudotime"] = X_plot
    
    step_number = 3
    for gene in list_of_genes:
        #Y_train = adata.obs_vector(gene, layer="imputed")
        Y_train = get_obs_vector(session_ID, gene, layer="imputed")

        gam = LinearGAM(n_splines=5, spline_order=3)
        gam.gridsearch(X_train, Y_train, weights=weights, progress=False)
        #gam = ExpectileGAM(terms="s(0)", expectile=0.5).gridsearch(X_train, Y_train)
        #lam = gam.lam
        #gam_upper = ExpectileGAM(expectile=0.75, lam=lam).fit(X_train, Y_train)
        #gam_lower = ExpectileGAM(expectile=0.25, lam=lam).fit(X_train, Y_train)
        gene_trends[gene] = gam.predict(X_plot)
        #gene_trends[gene + "_ci_upper"] = gam_upper.predict(X_plot)
        #gene_trends[gene + "_ci_lower"] = gam_lower.predict(X_plot)
        
        ci = gam.confidence_intervals(X_plot, width=.95)
        gene_trends[gene + "_ci_upper"] = ci[:,1]
        gene_trends[gene + "_ci_lower"] = ci[:,0]
        cache_progress(session_ID, progress=int(step_number/n_steps * 100))
        step_number += 1
    gene_trends = gene_trends.clip(lower = 0)
    return gene_trends