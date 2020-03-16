import dash
from dash.dependencies import Input, Output, State

from helper_functions import *
from plotting.plotting_functions import *
from app import app

from . pseudotime_functions import *

@app.callback(
	[Output("pseudotime_calculation_UMAP_plot", "figure"),
     Output("pseudotime_calculation_status", "children")],
	[Input("refresh_pseudotime_button", "n_clicks"),
	 Input("pseudotime_calculation_UMAP_dropdown", "value")],
	[State("session-id", "children"),
	 State("pseudotime_calculation_UMAP_plot", "selectedData")]
)
def refresh_pseudotime_UMAP(pt_btn_clicks, pt_plot_type, session_ID,
							selected_cells):
    print("[STATUS] refreshing pseudotime UMAP plot")
    default_return = [dash.no_update, dash.no_update]
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return default_return
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if(button_id == "refresh_pseudotime_button"):
        print("[DEBUG] refresh_pseudotime_button clicked")
        if (pt_btn_clicks in [None, 0]):
            return default_return
        
        adata = cache_adata(session_ID)
        if (adata is None):
            return default_return

        # get the starter cell from the pseudotime plot and redo the pseudotime based on that
        # if there are multiple selected, take the last one
        if not (selected_cells is None):

            pt_selected_cell_IDs = get_cell_intersection(session_ID, adata, [selected_cells])
                
        if ((pt_selected_cell_IDs is None) or (len(pt_selected_cell_IDs) != 1)):
            print("[ERROR] please select exactly 1 cell to use as a pseudotime starter cell"
                + "\nUsing first cell in set")
        starter_cell_ID = pt_selected_cell_IDs.pop()
        
        adata = do_pseudotime(session_ID, adata, starter_cell_ID)

    # if it's a dropdown menu update - load adata
    elif(button_id == "pseudotime_calculation_UMAP_dropdown"):
        if (pt_plot_type in [0, "", None, []]):
            return default_return
        else:
            adata = cache_adata(session_ID)

    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return default_return
    
    if (adata is None):
        return default_return

    # update the plot
    if (pt_plot_type in [0, "", None, []]):
        return default_return
    print("[STATUS] updating plot by: " + str(pt_plot_type))
    adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))

    if ("pseudotime" in adata.obs):
        pt_calc_status = "pseudotime calculated"
    else:
        pt_calc_status = ""

    if not (selected_cells is None):
        cells_to_highlight = get_cell_intersection(session_ID, adata, [selected_cells])
    else:
        cells_to_highlight = []
    if (pt_plot_type == "leiden_n"):
        return plot_UMAP(adata, "leiden", cells_to_highlight), pt_calc_status
    elif (pt_plot_type == "pseudotime"):
        return plot_pseudotime_UMAP(adata, "pseudotime"), pt_calc_status
    elif (pt_plot_type == "differentiation potential"):
        return plot_pseudotime_UMAP(adata, "differentiation_potential"), pt_calc_status
    elif (pt_plot_type == "total_counts"):
        return plot_expression_UMAP(adata, "total_counts"), pt_calc_status
    elif (pt_plot_type == "log1p_total_counts"):
        return plot_expression_UMAP(adata, "log1p_total_counts"), pt_calc_status
    elif (pt_plot_type == "n_genes"):
        return plot_expression_UMAP(adata, "n_genes"), pt_calc_status
    elif ("pseudotime_branch_" in pt_plot_type):
        return plot_pseudotime_UMAP(adata, pt_plot_type), pt_calc_status
    else:
        return default_return


# Refreshes dropdown options for pseudotime dropdown here 
# AND on the annotations page
@app.callback(
    [Output("pseudotime_dropdown", "options"),
     Output("pseudotime_gene_branch_dropdown", "options"),
     Output("pseudotime_calculation_UMAP_dropdown", "options")],
    [Input("pseudotime_calculation_status", "children")],
    [State("session-id", "children"),
     State("pseudotime_calculation_UMAP_dropdown", "options")]
)
def update_pseudotime_dropdown(pt_calc_status, session_ID, pt_calc_dropdown_options):
    default_return = [dash.no_update, dash.no_update, dash.no_update]

    if (pt_calc_status in [0, "", None, []]):
        return default_return

    adata = cache_adata(session_ID)
    if (adata is None):
        return default_return

    pt_branch_names = [col for col in adata.obs.columns if 'pseudotime_branch_' in col]
    value_list = ["pseudotime", "differentiation_potential"] + pt_branch_names
    options_0 = [{"label": i, "value": i} for i in value_list]
    options_1 = [{"label": i, "value": int(i[-1])} for i in pt_branch_names]
    options_2 = pt_calc_dropdown_options
    for x in [{"label": i, "value": i} for i in pt_branch_names]:
        unique = False
        for y in options_2:
            if (x == y):
               unique = False
               break
            else:
                unique = True
        if (unique == True):
            options_2.append(x)  
    return options_0, options_1, options_2