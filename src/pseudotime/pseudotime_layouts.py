import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from . import pseudotime_components as cc

def pseudotime_layout():
	m = dbc.Tab(label='Pseudotime', children=[
			dbc.Row(children=[
				dbc.Col(children=[
					# Buttons to do the recalculations
					dbc.Button("Recalculate pseudotime", 
				    			id="refresh_pseudotime_button"),
				    html.P('''
				    	Be sure to select a starter cell in the plot 
				    	below before running this.\nThis will take a
				    	lonnnngggg time and might fail to converge; please
				    	be patient.
						'''),
			    ], width=3),
			]),
			dbc.Row(children=[
				dbc.Col(children=[
					# pseudotime UMAP plot
				    html.H3(children="Pseudotime UMAP plots"),
				    html.P('''
				    	Use this dropdown menu to select a starting cell 
				    	before calculating a pseudotime trajectory, then use it
				    	to observe said trajectory.
				    	'''),
				    cc.pseudotime_UMAP_dropdown(),
			    	dcc.Loading(children=[cc.plot_pseudotime_UMAP()])		
				], width=6),
			]),
		]) # end pseudotime tab
	return m