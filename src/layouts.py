import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_dangerously_set_inner_html

from flask_security import current_user

import uuid
import base64

from processing.processing_layouts import *
from markergenes.markergenes_layouts import *
from pseudotime.pseudotime_layouts import *
from annotation.annotation_layouts import *
from importing.importing_layouts import *
from exporting.exporting_layouts import *
from status.status_layouts import *

demo = False 

def titlebar_layout(demo=False):	
	# load the logo
	MiCV_logo_filename = "/srv/www/MiCV/assets/MiCV_logo.png"
	MiCV_logo = base64.b64encode(open(MiCV_logo_filename, 'rb').read())

	if (demo == True):
		login_logout_link = ""
	elif (current_user and current_user.is_authenticated):
		login_logout_link = '<a class="nav-item nav-link" href="/logout">Logout</a>'
	else:
		login_logout_link = '<a class="nav-item nav-link" href="/login">Login</a>'

	ret = dash_dangerously_set_inner_html.DangerouslySetInnerHTML('''
		<nav class="navbar navbar-light bg-light">
          <a class="navbar-brand" href="/" className="ml-2">
            <img src="assets/MiCV_logo.png" height=45 alt="" loading="lazy">
            A Multi-informatic Cellular Visualization tool
          </a>
          <div class="navbar-expand flex-grow-1 text-left" id="navbarNav">
            <div class="navbar-nav"> ''' 
            + login_logout_link 
            + '''  
              <a class="nav-item nav-link" href="https://github.com/cailabumich/MiCV">Code</a>
              <a class="nav-item nav-link" href="https://micv.works/docs/">Documentation</a>
            </div>
          </div>
        </nav>
        ''')
	
	return ret

'''
	ret = html.Nav(children=[
			html.A(children=[
				html.Img(src="data:image/png;base64,{}".format(MiCV_logo.decode()), height="45px"),
				html.H3("A Multi-informatic Cellular Visualization tool"),
			], href="https://micv.works", className="navbar-brand"),
			html.Div(children=[
				dbc.Nav(children=[
					(html.A("Logout", href="/logout", className="nav-item nav-link") if (demo is False) else html.A("Login", href="/login", className="nav-item nav-link")),
					html.A("About", href="/about", className="nav-item nav-link"),
					html.A("Documentation", href="/documentation", className="nav-item nav-link")
			    ]),
			], className="navbar-expand flex-grow-1 text-left", id="navbarNav")
	], className="navbar navbar-light bg-light")

'''

def main_layout():
	if ((current_user) and (current_user.is_authenticated is True)):
		session_id = str(current_user.id)
	else:
		session_id = str(uuid.uuid4())

	return html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_1', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_2', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI

	    dcc.Interval(id="cleanup_interval",
	    			 interval=(15 * 60 * 1000),
	    			 max_intervals=-1),

		dcc.Interval(id="status_interval", 
					 interval=(1 * 3 * 1000),
					 max_intervals=-1),
		# title matter
		titlebar_layout(demo=demo),

	    dbc.Container(fluid=True, children=[
		    dbc.Row(children=[
		    	# Main layout
		    	dbc.Col(children=[
		    		dbc.Tabs(children=[
		    			### IMPORTING TAB ###
			        	dbc.Tab(importing_layout(demo=demo), label="Load data", tab_id="importing_tab"),

				    	### PROCESSING TAB ###
			        	dbc.Tab(processing_layout(demo=demo), label="Preprocess", tab_id="processing_tab"),
						
						### MARKER GENE TAB ###
						dbc.Tab(markergenes_layout(), label="Marker genes", tab_id="markergenes_tab"),

						### PSEUDOTIME TAB ###
						dbc.Tab(pseudotime_layout(), label="Pseudotime", tab_id="pseudotime_tab"),
						
					    ### ANNOTATION TAB ###
					    dbc.Tab(annotation_layout(), label="Expert annotation", tab_id="annotation_tab"),

						### DOWNLOAD ANALYSIS TAB ###
						dbc.Tab(exporting_layout(demo=demo), label="Save and export", tab_id="exporting_tab"),
					], id="main_tabs"),	
		    	], width=10),

		    	# Status panel
		    	dbc.Col(children=[
		    		status_layout()
		    	], width=2)
		    ])

		]),
	]) # end main layout
