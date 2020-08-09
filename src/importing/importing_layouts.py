import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from . import importing_components as cc

def importing_layout(demo=False):
	m = dbc.Container(children=[
			dbc.Row(children=[
				dbc.Col(children=[
					cc.importing_greeting(demo=demo)
				])
			]),
			dbc.Row(children=[
		    	dbc.Col(children=[
					dbc.CardDeck(children=[
						cc.importing_dataset_dropdown(),
						cc.importing_user_dataset_list(demo=demo),
			        	cc.importing_data_upload(demo=demo),
		        	]),
			    ])
			], id="upload-collapse")
		], fluid=True)
	return m