import scanpy as sc

def identify_marker_genes(adata, obs_column, groups_to_rank, method):
    if (method == "logreg"):
	    if ("all" in groups_to_rank):
		    if (len(adata.obs[obs_column].unique()) < 3):
			    print("[ERROR] too few groups to use logreg maker gene detection;"+
					  " using t-test_overestim_var by default")
			    method = "t-test_overestim_var"
	    else:
		    if (len(groups_to_rank) < 3):
			    print("[ERROR] too few groups to use logreg maker gene detection;"+
					  " using t-test_overestim_var by default")
			    method = "t-test_overestim_var"

    if ("all" in groups_to_rank):
        sc.tl.rank_genes_groups(adata, groupby=obs_column, method=method)
    else:
        sc.tl.rank_genes_groups(adata, groupby=obs_column, 
                                groups=groups_to_rank, method=method)
    return adata