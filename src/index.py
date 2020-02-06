import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import uuid

from app import app
import layouts
import callbacks

#session_id = str(uuid.uuid4())
session_id = "testing"

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/MiCV':
        session_id, layout =  layouts.main_layout()
        return layout
    elif pathname == '/information':
        session_id, layout = layouts.information()
        return layout
    else:
        return "404 URL not found"

server = app.server

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')