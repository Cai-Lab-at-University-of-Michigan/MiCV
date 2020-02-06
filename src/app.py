import dash
from flask_caching import Cache
import layouts


app = dash.Dash(__name__, show_undo_redo=False)

app.server.secret_key = b"\xc9j\xa2@k\x04\x0e\x8a\xe9\xb6\xfbA\xdfsU\x05\xdfe\xec@\x05\x0b\xfd\x9a"
server = app.server

app.config.suppress_callback_exceptions = True

'''
cache = Cache(app.server, config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': '/tmp/MiCV',
    'CACHE_THRESHOLD': 200,  # should be equal to maximum number of active users
    "CACHE_DEFAULT_TIMEOUT": 30000
    }
)
'''