import dash_html_components as html
import dash_bootstrap_components as dbc

from . import markergenes_components as cc #custom components

def markergenes_layout():
	m = dbc.Tab(label="Marker genes", children=[
				dbc.Row(children=[
					dbc.Col(children=[
						html.H3(children="Clustering plot"),
						html.P(
							'''
							obs columns 
							can be used to define groups for marker gene 
							identification (typically you should use the 
							leiden column or one of the user_n manually 
							annotated columns)
							'''),
						cc.marker_gene_UMAP_dropdown(),
						cc.plot_marker_gene_UMAP()
					], width=6),
				]),
				dbc.Row(children=[
					dbc.Col(children=[
						html.H3(children="Marker gene identification"),
						html.P(
							'''
							Select which groups you would like to identify
							marker genes for and which method you would like
							to use for identification
							'''),
						cc.marker_gene_group_dropdown(),
						cc.marker_gene_method_radio(),
						dbc.Button("Recalculate marker genes", 
					    			id="recalculate_marker_genes"),
						cc.marker_gene_plot()
					]),
				])
			]) # end marker gene tab
	return m
