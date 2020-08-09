import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

from . import markergenes_components as cc #custom components

def markergenes_layout():
	m = dbc.Tab(label="Marker genes", children=[
				dbc.Row(children=[
					dbc.Col(children=[
						html.H3(children="Clustering plot"),
						html.P(
							'''
							Select a categorical observation (cell)
							column that contains the categories you
							would like to calculate marker genes 
							for. Then, in the dropdown menu below
							this plot, select which categories you 
							would like to compare.
							'''),
						cc.marker_gene_UMAP_dropdown(),
						cc.n_dims_proj_markergenes_radio(),
						dcc.Loading(children=[cc.plot_marker_gene_UMAP()])
					]),
				]),
				dbc.Row(children=[
					dbc.Col(children=[
						html.H3(children="Marker gene identification"),
						html.P(
							'''
							Select the groups you would like to compare,
							and select which marker gene detection algorithm 
							you would like to use for said comparison. 
							Logistic-regression gives 
							very meaningful results, but takes longer
							to run and requires that you compare at
							least 3 groups at a time. Consider
							t-test_overestim_var for a quick, but 
							statistically meaningful, look at 
							putative marker genes.
							'''),
						cc.marker_gene_group_dropdown(),
						cc.marker_gene_method_radio(),
						dbc.Row(children=[
							dbc.Col(children=[
								dbc.Button("Recalculate marker genes", 
					    					id="recalculate_marker_genes")
							], width=2),
							dbc.Col(children=[
								cc.marker_genes_export_button()
							], width=2)
						]),
						dcc.Loading(children=[cc.marker_gene_plot()])
					], width=12)
				])
			]) # end marker gene tab
	return m
