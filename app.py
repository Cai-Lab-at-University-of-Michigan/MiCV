import dash
from flask_caching import Cache
import layouts

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
				show_undo_redo=False)
app.layout = layouts.main_layout
server = app.server

app.config.suppress_callback_exceptions = True

cache = Cache(app.server, config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': '/tmp/MiCV',
    'CACHE_THRESHOLD': 200,  # should be equal to maximum number of active users
    "CACHE_DEFAULT_TIMEOUT": 30000
    }
)
