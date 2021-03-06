import os
import shutil
import time
import pickle
from datetime import datetime

from filelock import Timeout, FileLock
import zarr
from numcodecs import Blosc

import pandas as pd
import numpy as np
import scipy as sp
import scanpy as sc
import anndata as ad
from anndata._io.zarr import read_dataframe, read_attribute, write_attribute

from status.status_functions import *
from plotting.multi_color_scale import MultiColorScale

from tasks.tasks import write_dense

save_analysis_path = "/srv/www/MiCV/cache/"
selected_datasets_path = "/srv/www/MiCV/selected_datasets/"
user_dataset_path  = "/srv/www/MiCV/user_datasets/"

lock_timeout = 60

use_zarr = True

### the actual helper functions
### TODO: break all of this up into specific modules

def generate_adata_from_10X(session_ID, data_type="10X_mtx"):
    data_dir = save_analysis_path + str(session_ID) + "/raw_data/"
    if (data_type == "10X_mtx"):
        adata = sc.read_10x_mtx(data_dir, cache=False)
    elif (data_type == "10X_h5"):
        adata = sc.read_10x_h5(data_dir + "data.h5ad")
    else:
        print("[ERROR] data type not recognized - returning None")
        return None

    cache_adata(session_ID, adata)
    return adata

def load_selected_dataset(session_ID, dataset_key):
    dataset_dict = {
    "00001": "Michki2020",
    "00002": "Cocanougher2020",
    "00003": "Davie2018",
    "00004": "10X5KPBMC",
    "00005": "Sharma2020",
    "00006": "Zeisel2018"
    }

    filename = dataset_dict[dataset_key]
    if (filename is None):
        return None
    else:
        filename = selected_datasets_path + filename

    adata = sc.read_h5ad(filename + ".h5ad")

    state = {"filename": str(dataset_dict[dataset_key]),
             "# cells/obs": len(adata.obs.index),
             "# genes/var": len(adata.var.index),
             "# counts": int(np.sum(adata.obs["total_counts"]))}
    cache_state(session_ID, state)

    adata = cache_adata(session_ID, adata)
    return adata

def cache_adata(session_ID, adata=None, group=None,
                store_dir=None, store_name=None):
    if ((store_dir is None) or (store_name is None)):
        save_dir = save_analysis_path  + str(session_ID) + "/"
        filename = save_dir + "adata_cache"
        chunk_factors = [150, 3] #faster, hot storage
    else:
        save_dir = store_dir
        filename = save_dir + store_name
        chunk_factors = [3, 3] #slower, cold storage

    
    if not (os.path.isdir(save_dir)):
        try:
            print("[DEBUG] making directory:" + str(save_dir))
            os.mkdir(save_dir)
        except:
            return None
    
    lock_filename = (save_analysis_path  + str(session_ID) 
                     + "/" + "adata.lock")
    lock = FileLock(lock_filename, timeout=lock_timeout)

    compressor = Blosc(cname='blosclz', clevel=3, 
                       shuffle=Blosc.SHUFFLE)
    zarr_cache_dir = filename  + ".zarr"
    attribute_groups = ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "X", "uns", "raw"]
    extra_attribute_groups = ["X_dense", "layers_dense"]

    if (adata is None): # then -> read it from the store
        if (os.path.exists(zarr_cache_dir) is True):
            store_store = zarr.DirectoryStore(zarr_cache_dir)
            store = zarr.open_group(store=store_store, mode='r')
            if (group in attribute_groups): # then -> return only that part of the object (fast)
                group_exists = adata_cache_group_exists(session_ID, group, store=store)
                if (group_exists is True):
                    ret = read_attribute(store[group])
                else:
                    ret = None
                #store_store.close()
                return ret            
            elif (group is None): # then -> return the whole adata object (slow)
                #adata = ad.read_zarr(zarr_cache_dir)
                d = {}
                for g in attribute_groups:
                    if (g in store.keys()):
                        if (adata_cache_group_exists(session_ID, g, store=store)):
                            if (g in ["obs", "var"]):
                                d[g] = read_dataframe(store[g])
                            else:
                                d[g] = read_attribute(store[g])
                #store_store.close()
                adata = ad.AnnData(**d)
                if not (adata is None):
                    return adata
                else:
                    print("[ERROR] adata object not saved at: " + str(filename))
                    return None
    else: # then -> update the state dictionary and write adata to the store
        if (group is None):
            cache_state(session_ID, key="# cells/obs", val=len(adata.obs.index))
            cache_state(session_ID, key="# genes/var", val=len(adata.var.index))
            if ("total_counts" in adata.obs):
                cache_state(session_ID, key="# counts", val=int(np.sum(adata.obs["total_counts"])))
            else:
                cache_state(session_ID, key="# counts", val=int(np.sum(adata.X)))

        elif (group == "obs"):
            cache_state(session_ID, key="# cells/obs", val=len(adata.index))
        elif (group == "var"):
            cache_state(session_ID, key="# genes/var", val=len(adata.index))
        with lock:
            store_store = zarr.DirectoryStore(zarr_cache_dir)
            store = zarr.open_group(store=store_store, mode='a')
            if (group in attribute_groups): # then -> write only that part of the object (fast)
                if (group == "var"):
                    if (np.nan in adata.var.index):
                        adata.var.index = pd.Series(adata.var.index).replace(np.nan, 'nanchung')
                        adata.var["gene_ID"] = pd.Series(adata.var["gene_ID"]).replace(np.nan, 'nanchung')
                        adata.var["gene_ids"] = pd.Series(adata.var["gene_ids"]).replace(np.nan, 'nanchung')
                write_attribute(store, group, adata) # here "adata" is actually just a subset of adata
                
                # write dense copies of X or layers if they're what was passed
                if (group == "X"):
                    dense_name = "X_dense"
                    write_dense.delay(zarr_cache_dir, "X",
                                      dense_name, chunk_factors)

                if (group == "layers"):
                    for l in list(adata.keys()): #layers was passed with parameter name "adata"
                        dense_name = "layers_dense/" + str(l)
                        write_dense.delay(zarr_cache_dir, "layers/" + l, 
                                          dense_name, chunk_factors)
                #store_store.flush()
                #store_store.close()
                lock.release()
            else:
                # check that necessary fields are present in adata object
                if not ("leiden_n" in adata.obs):
                    if ("leiden" in adata.obs):
                        adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])
                if not ("cell_ID" in adata.obs):
                    adata.obs["cell_ID"] = adata.obs.index
                if not ("cell_numeric_index" in adata.obs):
                    adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
                for i in ["user_" + str(j) for j in range(0, 6)]:
                    if not (i in adata.obs.columns):
                        adata.obs[i] = ["0" for k in adata.obs.index.to_list()]
                if not ("gene_ID" in adata.var):
                    adata.var["gene_ID"] = adata.var.index

                # make sure that there are no "nan" genes in the var index
                if (np.nan in adata.var.index):
                    adata.var.index = pd.Series(adata.var.index).replace(np.nan, 'nanchung')
                    adata.var["gene_ID"] = pd.Series(adata.var["gene_ID"]).replace(np.nan, 'nanchung')
                    adata.var["gene_ids"] = pd.Series(adata.var["gene_ids"]).replace(np.nan, 'nanchung')

                # save it all to the cache, but make dense copies of X and layers
                write_attribute(store, "obs", adata.obs)
                write_attribute(store, "var", adata.var)
                write_attribute(store, "obsm", adata.obsm)
                write_attribute(store, "varm", adata.varm)
                write_attribute(store, "obsp", adata.obsp)
                write_attribute(store, "varp", adata.varp)
                write_attribute(store, "uns", adata.uns)
                write_attribute(store, "raw", adata.raw)
                write_attribute(store, "X", adata.X)
                write_attribute(store, "layers", adata.layers)

                # making dense copies of X and layers (compressed to save disk space)
                
                dense_name = "X_dense"
                write_dense.delay(zarr_cache_dir, "X",
                                  dense_name, chunk_factors)

                for l in list(adata.layers.keys()):
                    dense_name = "layers_dense/" + str(l)
                    write_dense.delay(zarr_cache_dir, "layers/" + l, 
                                      dense_name, chunk_factors)
                
                #store_store.flush()
                #store_store.close()
                lock.release()
            # set the file mod and access times to current time
            # then return adata as usual 
            os.utime(zarr_cache_dir) 
            return adata

def adata_cache_exists(session_ID):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    filename = save_dir + "adata_cache"
    zarr_cache_dir = filename  + ".zarr"
        
    if (os.path.exists(zarr_cache_dir) is True):
        return True
        
    return False

def adata_cache_group_exists(session_ID, group, store=None):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    filename = save_dir + "adata_cache"
    zarr_cache_dir = filename  + ".zarr"
    
    if (store is None):
        try:
            store_store = zarr.DirectoryStore(zarr_cache_dir)
            store = zarr.open_group(store=store_store, mode='r')
            keys = list(store.group_keys())
            #store_store.close()
        except:
            return False
    else:
        keys = list(store.group_keys())
    
    if (group in keys):
        return True
    else:
        return False

def cache_gene_list(session_ID, gene_list=None):
    filename = save_analysis_path + str(session_ID) + "/gene_list_cache.pickle"
    lock_filename = filename + ".lock"
    lock = FileLock(lock_filename, timeout=20)
    
    if (gene_list is None):
        if (os.path.isfile(filename) is True):
            with open(filename, "rb") as f:
                gene_list = pickle.load(f)
        else:
            print("[ERROR] gene list cache does not exist at: " + str(filename))
            gene_list = None
        return gene_list
    else:
        gene_list.sort(key=str.lower)
        with lock:
            with open(filename, "wb") as f:
                pickle.dump(gene_list, f)
            return gene_list

# returns a list of cell_IDs 
# expects a list of lists of datapoint dictionaries
def get_cell_intersection(session_ID, list_of_selections,
                          pt_min=0, pt_max=1):
    obs = cache_adata(session_ID, group="obs")
    cell_intersection = set(obs.index.to_list())

    for cell_list in list_of_selections:
        if (cell_list in ["", 0, None, []]):
            continue
        
        cell_set = set()
        for cell in cell_list["points"]:
            cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
            cell_set.add(cell_ID)
        cell_intersection &= cell_set
    
    if ("pseudotime" in obs): 
        if ((pt_min > 0) or (pt_max < 1)):
            for cell in list(cell_intersection):
                if ((obs.loc[cell, "pseudotime"] < pt_min)
                or  (obs.loc[cell, "pseudotime"] > pt_max)):
                    cell_intersection.remove(cell)

    return cell_intersection

# returns dictionary of points that are in all violin selections
# across all genes that could be selected on
def get_violin_intersection(session_ID, violin_selected):
    if (violin_selected is None):
        return None

    # test which traces (genes) cells were selected from
    curves = set()
    for cell in violin_selected["points"]:
        curves.add(cell["curveNumber"])

    if (len(curves) == 0):
        # no cells selected - return None
        return None

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

def cache_multicolor_scale(session_ID, multi_color_scale=None):
    if not (session_ID == None):
        filename = save_analysis_path + str(session_ID) + "/multi_color_scale.pickle"
    else:
        filename = save_analysis_path + "multi_color_scale.pickle"

    if (multi_color_scale is None):
        if (os.path.isfile(filename) is True):
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

def remove_old_cache(n_days=1.5):
    n_sec_in_day = 86400
    max_time_in_sec = int(n_sec_in_day * n_days)
    
    now = time.time()

    for r,d,f in os.walk(save_analysis_path):
        for directory in d:
            timestamp = os.path.getmtime(os.path.join(r,directory))
            if (now-max_time_in_sec > timestamp):
                try:
                    print("[CLEANUP] removing " + str(os.path.join(r,directory)))
                    shutil.rmtree(os.path.join(r,directory))
                except Exception as e:
                    print("[ERROR] " + str(e))
                    pass

# implementation specific to custom zarr file cache with X_dense and layers_dense
def get_obs_vector(session_ID, var, layer="X"):
    save_dir = save_analysis_path  + str(session_ID) + "/"

    if (use_zarr is True):
        zarr_cache_dir = save_dir + "adata_cache" + ".zarr"
        if (os.path.exists(zarr_cache_dir) is True):
            with zarr.LMDBStore(zarr_cache_dir) as store_store:
                store = zarr.open_group(store=store_store, mode='r')

                if (var in store.obs.keys()):
                    ret = list(store.obs[var])
                else:
                    idx = list(store.var["gene_ID"]).index(var)
                    if (layer == "X"):
                        ret = store["X_dense"][:, idx]
                    else:
                        ret = (store["layers_dense"][layer])[:, idx]
            return ret

def generate_marker_gene_table(session_ID):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    uns = cache_adata(session_ID, group="uns")
    marker_genes = pd.DataFrame.from_records(uns["rank_genes_groups"]["names"])
    marker_genes.to_csv(save_dir + "marker_genes.csv")

def marker_genes_table_exists(session_ID):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    filename = save_dir + "marker_genes.csv"
    if (os.path.isfile(filename) is True):
        return True
    return False

def to_rgba_string(rgb_tuple, opacity=1):
    ret = "rgba("
    for c in rgb_tuple:
        ret += str(c) + ","
    ret += str(opacity) + ")"
    return ret