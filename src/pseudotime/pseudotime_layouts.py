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
					html.H3(children="Pseudotime trajectory calculation"),
				    html.P('''
				    	Select a starter cell in the plot 
				    	below before running this. Expect it to take between 
				    	3 and 30 minutes to run, and recognize that
				    	it might fail to converge altogether. Be prepared 
				    	to re-run it with a different starter cell if 
				    	it fails to converge. 

				    	Consider skipping this tool altogether during
				    	preliminary data analysis.
						'''),
					dbc.Button("Recalculate pseudotime", 
				    			id="refresh_pseudotime_button"),
				    html.Div(
				    	 id="pseudotime_calculation_status",
				    	 children=[]
				    )
			    ], width=4),
			]),
			dbc.Row(children=[
				dbc.Col(children=[
					# pseudotime UMAP plot
				    html.H3(children="Pseudotime UMAP plots"),
				    html.P('''
				    	Use this dropdown menu to select a starting cell 
				    	before calculating a pseudotime trajectory, then use 
				    	it to observe said trajectory.
				    	'''),
				    cc.pseudotime_UMAP_dropdown(),
			    	dcc.Loading(children=[cc.plot_pseudotime_UMAP()])		
				]),
			]),
		]) # end pseudotime tab
	return m