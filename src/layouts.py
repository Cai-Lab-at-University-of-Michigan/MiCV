import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import uuid
import base64

from processing.processing_layouts import *
from markergenes.markergenes_layouts import *
from pseudotime.pseudotime_layouts import *
from annotation.annotation_layouts import *
from exporting.exporting_layouts import *
from status.status_layouts import *

def main_layout():
	session_id = str(uuid.uuid4())

	# load the logo
	MiCV_logo_filename = "/srv/www/MiCV/assets/MiCV_logo.png"
	MiCV_logo = base64.b64encode(open(MiCV_logo_filename, 'rb').read())

	return html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_1', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_2', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI

	    dcc.Interval(id="cleanup_interval",
	    			 interval=(15 * 60 * 1000),
	    			 max_intervals=-1),

		dcc.Interval(id="status_interval", 
					 interval=(1 * 10 * 1000),
					 max_intervals=-1),

	    dbc.Container(fluid=True, children=[
		    # title matter
		    dbc.Navbar(children=[
		    	html.A(
		    		dbc.Row(children=[
		    			dbc.Col(html.Img(src="data:image/png;base64,{}".format(MiCV_logo.decode()), height="45px")),
						dbc.Col(dbc.NavbarBrand("A Multi-informatic Cellular Visualization tool", className="ml-2"))
		    		], align="center", no_gutters=True),
		    	href="https://micv.works")
		    ]),

		    dbc.Row(children=[
		    	# Main layout
		    	dbc.Col(children=[
		    		dbc.Tabs(children=[
				    	### PROCESSING TAB ###
			        	dbc.Tab(processing_layout(), label="Load and preprocess", tab_id="processing_tab"),
						
						### MARKER GENE TAB ###
						dbc.Tab(markergenes_layout(), label="Marker genes", tab_id="markergenes_tab"),

						### PSEUDOTIME TAB ###
						dbc.Tab(pseudotime_layout(), label="Pseudotime", tab_id="pseudotime_tab"),
						
					    ### ANNOTATION TAB ###
					    dbc.Tab(annotation_layout(), label="Expert annotation", tab_id="annotation_tab"),

						### DOWNLOAD ANALYSIS TAB ###
						dbc.Tab(exporting_layout(), label="Save and export", tab_id="exporting_tab"),
					], id="main_tabs"),	
		    	], width=10),

		    	# Status panel
		    	dbc.Col(children=[
		    		status_layout()
		    	], width=2)
		    ])

		]),
]) # end main layout