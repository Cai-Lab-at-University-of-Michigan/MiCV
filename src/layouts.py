import dash_html_components as html
import dash_bootstrap_components as dbc
import uuid

from processing.processing_layouts import *
from markergenes.markergenes_layouts import *
from annotation.annotation_layouts import *
from exporting.exporting_layouts import *

def main_layout(session_id = None):
	if (session_id is None):
		session_id = "testing"
		#session_id = str(uuid.uuid4())

	return html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_1', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI

	    dbc.Container(fluid=True, children=[
		    # title matter
		    html.H1(children="MiCV"),
		    html.H3(children="A Multi-Informatic Cellular Visualization tool"),

		    dbc.Tabs(children=[
		    	### PROCESSING TAB ###
	        	processing_layout(),
				
				### MARKER GENE TAB ###
				markergenes_layout(),
				
			    ### ANNOTATION TAB ###
			    annotation_layout(),

				### DOWNLOAD ANALYSIS TAB ###
				exporting_layout(),
			]),
		]),
]) # end main layout