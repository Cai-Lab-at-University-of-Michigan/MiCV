import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
import layouts
import callbacks

app.layout = layouts.main_layout
server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)