import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import flask

import uuid

from app import app
import layouts

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

from inputoutput import inputoutput_callbacks
from processing import processing_callbacks
from markergenes import markergenes_callbacks
from pseudotime import pseudotime_callbacks
from annotation import annotation_callbacks
from exporting import exporting_callbacks
from status import status_callbacks

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/':
        layout =  layouts.main_layout()
        return layout
    else:
        return "404 URL not found"

server = app.server

# extra static routes
@server.route('/favicon.ico')
def send_favicon():
    return flask.send_from_directory("/srv/www/MiCV/assets/", 
    								 "favicon.ico")


if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')