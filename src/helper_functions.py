from os import path
import pickle
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np

from plotting.multi_color_scale import MultiColorScale

save_analysis_path = "/srv/www/MiCV/cache/"


def generate_adata_from_10X(session_ID, data_type="10X_mtx"):
    #data_dir = "/home/nigelmic/bioinformatics/Solo.out"
    data_dir = save_analysis_path + "/" + "raw_data"
    print("[STATUS] loading data from " + str(data_dir))
    if (data_type == "10X_mtx"):
        adata = sc.read_10x_mtx(data_dir, cache=True)
    elif (data_type == "10X_h5"):
        adata = sc.read_10x_h5(data_dir + "data.h5ad")
    else:
        print("[ERROR] data type not recognized - returning None")
        return None


    cache_adata(session_ID, adata)
    return adata

def cache_adata(session_ID, adata=None):
    #filename = save_analysis_path  + "adata_cache.h5ad"
    filename = save_analysis_path  + "adata_cache"
    #print("[DEBUG] filename = " + str(filename))
    if (adata is None):
        if (path.isfile(filename  + ".h5ad") is True):
            adata = sc.read_h5ad(filename + ".h5ad")

        if not (adata is None):
            if not ("leiden_n" in adata.obs):
                adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])
            if not ("cell_ID" in adata.obs):
                adata.obs["cell_ID"] = adata.obs.index
            if not ("cell_ID" in adata.obs):
                adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
            adata.write(filename + ".h5ad")
            #adata.write_zarr(filename + ".zarr", (1000,1000))
            adata.filename = None
            #return ad.AnnData(X = adata.X, var=adata.var, obs=adata.obs, uns=adata.uns, obsm=adata.obsm, varm=adata.varm, raw=adata.raw)
            return adata
        else:
            print("[ERROR] adata object not saved at: " + str(filename))
            return None
    else:
        adata.write(filename + ".h5ad")
        #adata.write_zarr(filename + ".zarr", (1000,1000))
        adata.filename = None
        #return ad.AnnData(X = adata.X, var=adata.var, obs=adata.obs, uns=adata.uns, obsm=adata.obsm, varm=adata.varm, raw=adata.raw)
        return adata

def cache_gene_trends(session_ID, gene_trends=None):
    filename = save_analysis_path + "gene_trends_cache.pickle"
    if (gene_trends is None):
        if (path.isfile(filename) is True):
            with open(filename, "rb") as f:
                gene_trends = pickle.load(f)
        else:
            print("[ERROR] gene trends cache does not exist at: " + str(filename))
            gene_trends = ["NULL"]
        return gene_trends
    else:
        with open(filename, "wb") as f:
            pickle.dump(gene_trends, f)
        return gene_trends

def cache_gene_list(session_ID, gene_list=None):
    filename = save_analysis_path + "gene_list_cache.pickle"
    if (gene_list is None):
        if (path.isfile(filename) is True):
            with open(filename, "rb") as f:
                gene_list = pickle.load(f)
        else:
            print("[ERROR] gene list cache does not exist at: " + str(filename))
            gene_list = ["NULL"]
        return gene_list
    else:
        with open(filename, "wb") as f:
            pickle.dump(gene_list, f)
        return gene_list

# returns a list of cell_IDs 
# expects a list of lists of datapoint dictionaries
def get_cell_intersection(session_ID, adata, list_of_selections,
                          pt_min=0, pt_max=1):
    cell_intersection = set(adata.obs.index.to_list())

    for cell_list in list_of_selections:
        if (cell_list in ["", 0, None, []]):
            continue
        
        cell_set = set()
        for cell in cell_list["points"]:
            cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
            cell_set.add(cell_ID)
        cell_intersection &= cell_set
    
    if ("pseudotime" in adata.obs): 
        if ((pt_min > 0) or (pt_max < 1)):
            for cell in list(cell_intersection):
                if ((adata.obs.loc[cell, "pseudotime"] < pt_min)
                or  (adata.obs.loc[cell, "pseudotime"] > pt_max)):
                    cell_intersection.remove(cell)

    return cell_intersection

# returns dictionary of points that are in all violin selections
# across all genes that could be selected on
def get_violin_intersection(session_ID, adata, violin_selected):
    if (violin_selected is None):
        return None

    # test which traces (genes) cells were selected from
    curves = set()
    for cell in violin_selected["points"]:
        curves.add(cell["curveNumber"])

    # get cells selected in each of these curves
    cells_in_curves = [set() for curve in curves]
    
    if (len(cells_in_curves) == 1):
        return violin_selected

    for cell in violin_selected["points"]:
        n = cell["curveNumber"]
        if (n in curves):
            cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
            (cells_in_curves[n]).add(cell_ID)

    # do the intersection
    cell_intersection = cells_in_curves[0]
    for cell_set in cells_in_curves:
        cell_intersection &= cell_set 

    # get the final list of points in dict format
    points = []
    for cell in violin_selected["points"]:
        cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
        if (cell_ID in cell_intersection):
            points.append(cell)

    return {"points": points}

def get_pseudotime_min_max(session_ID, selected):
    if (selected in ["", 0, None, []]):
            return [0,1]        
    
    x_vals = []
    for point in selected["points"]:
        x_vals.append(point["x"])

    if (len(x_vals) == 0):
        return [0,1]

    x_min = np.min(x_vals)
    x_max = np.max(x_vals)
    
    return [x_min, x_max]

def get_ortholog_data(session_ID, selected_gene):
    filename = (save_analysis_path 
    + "dmel_human_orthologs_disease_fb_2019_05.csv")

    ortholog_data = pd.read_csv(filename, sep="\t")
    ret = ortholog_data.loc[ortholog_data["Dmel_gene_symbol"] == selected_gene]
    return ret

def get_gene_snapshot(session_ID, selected_gene):
    filename = (save_analysis_path 
    + "gene_snapshots_fb_2019_05.csv")

    snapshots = pd.read_csv(filename, sep="\t")
    ret = snapshots.loc[snapshots["GeneSymbol"] == selected_gene]
    #ret = ret["gene_snapshot_text"]
    return ret

def get_disease_data(session_ID, selected_gene):
    filename = (save_analysis_path 
    + "disease_model_annotations_fb_2019_05.csv")

    diseases = pd.read_csv(filename, sep="\t")
    ret = diseases.loc[diseases["Dmel_gene_ID"] == selected_gene]
    #ret = ret["gene_snapshot_text"]
    return ret

def cache_multicolor_scale(multi_color_scale=None):
    filename = save_analysis_path + "multi_color_scale.pickle"
    if (multi_color_scale is None):
        if (path.isfile(filename) is True):
            with open(filename, "rb") as f:
                multi_color_scale = pickle.load(f)
        else:
            print("[ERROR] multi_color_scale cache does not exist at: " + str(filename))
            multi_color_scale = MultiColorScale()
            with open(filename, "wb") as f:
                pickle.dump(multi_color_scale, f)
        return multi_color_scale
    else:
        with open(filename, "wb") as f:
            pickle.dump(multi_color_scale, f)
        return multi_color_scale