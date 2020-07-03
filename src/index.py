import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_dangerously_set_inner_html
from dash.dependencies import Input, Output

import flask
import flask_security
from flask_security import current_user, logout_user

import uuid

from app import app
import layouts

app.layout = layouts.main_layout

from inputoutput import inputoutput_callbacks
from processing import processing_callbacks
from markergenes import markergenes_callbacks
from pseudotime import pseudotime_callbacks
from annotation import annotation_callbacks
from importing import importing_callbacks
from exporting import exporting_callbacks
from status import status_callbacks

server = app.server

# extra static routes
@server.route('/favicon.ico')
def send_favicon():
    return flask.send_from_directory("/srv/www/MiCV/assets/", 
    								 "favicon.ico")

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')