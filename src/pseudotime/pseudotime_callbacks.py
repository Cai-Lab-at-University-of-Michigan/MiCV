import dash
from dash.dependencies import Input, Output, State

from helper_functions import *
from plotting_functions import *
from app import app

from . pseudotime_functions import *

@app.callback(
	Output("pseudotime_calculation_UMAP_plot", "figure"),
	[Input("refresh_pseudotime_button", "n_clicks"),
	 Input("pseudotime_calculation_UMAP_dropdown", "value")],
	[State("session-id", "children"),
	 State("pseudotime_calculation_UMAP_plot", "selectedData")]
)
def refresh_pseudotime_UMAP(pt_btn_clicks, pt_plot_type, session_ID,
							selected_cells):
    print("[STATUS] refreshing pseudotime UMAP plot")
    # figure out which button was pressed - what refresh functions to call
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "not_triggered"
        return dash.no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if(button_id == "refresh_pseudotime_button"):
        print("[DEBUG] refresh_pseudotime_button clicked")
        if (pt_btn_clicks in [None, 0]):
            return dash.no_update
        
        # get the starter cell from the pseudotime plot and redo the pseudotime based on that
        # if there are multiple selected, take the last one
        if not (selected_cells is None):
            adata = cache_adata(session_ID)
            pt_selected_cell_IDs = get_cell_intersection(session_ID, adata, [selected_cells])
        
        if ((pt_selected_cell_IDs is None) or (len(pt_selected_cell_IDs) != 1)):
            print("[ERROR] please select exactly 1 cell to use as a pseudotime starter cell"
                + "\nUsing first cell in set")
        starter_cell_ID = pt_selected_cell_IDs.pop()
        
        adata = do_pseudotime(session_ID, adata, starter_cell_ID)

    # if it's a dropdown menu update - load adata
    elif(button_id == "pseudotime_calculation_UMAP_dropdown"):
        if (pt_plot_type in [0, "", None, []]):
            return dash.no_update
        else:
            adata = cache_adata(session_ID)


    # do nothing if no buttons pressed
    elif(button_id == "not_triggered"):
        return dash.no_update
    
    # update the plot
    if (pt_plot_type in [0, "", None, []]):
        return dash.no_update
    print("[STATUS] updating plot by: " + str(pt_plot_type))
    adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))

    if not (selected_cells is None):
        cells_to_highlight = get_cell_intersection(session_ID, adata, [selected_cells])
    else:
        cells_to_highlight = []
    if (pt_plot_type == "leiden_n"):
        return plot_UMAP(adata, "leiden_n", cells_to_highlight)
    elif (pt_plot_type == "pseudotime"):
        return plot_pseudotime_UMAP(adata, "pseudotime")
    elif (pt_plot_type == "differentiation potential"):
        return plot_pseudotime_UMAP(adata, "differentiation_potential")
    elif (pt_plot_type == "total_counts"):
        return plot_expression_UMAP(adata, "total_counts")
    elif (pt_plot_type == "log1p_total_counts"):
        return plot_expression_UMAP(adata, "log1p_total_counts")
    elif (pt_plot_type == "n_genes"):
        return plot_expression_UMAP(adata, "n_genes")