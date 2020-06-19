import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from . import importing_components as cc

def importing_layout(demo=False):
	m = dbc.Row(children=[
	    	dbc.Col(children=[
	    		html.Div([
				    html.P('''
				    	Either upload your own data (disabled) 
				    	or search for a pre-made dataset
				    	to load into MiCV for analysis.
						''')
					], style={'marginBottom': 20, 'marginTop': 20}),
		        dbc.Row(children=[
		        	cc.importing_dataset_dropdown(),
		        	(cc.importing_data_upload() if (demo is False) else html.Div()),
				], no_gutters=True),
				html.Div(id='upload_spacer_div',
			    		 style={'margin-top': 250}),			
			], id="upload-collapse")
		])

	return m