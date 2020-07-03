import dash

from flask import Flask
from flask_mail import Mail
from flask_caching import Cache
from flask_security import Security, login_required, \
     SQLAlchemySessionUserDatastore
from sqlalchemy import exc
from flask_bootstrap import Bootstrap

from user_management.database import db_session, init_db
from user_management.models import User, Role

import layouts
from tasks.tasks import send_flask_mail
from configmodule.default_config import DefaultConfig as FlaskConfig

security = Security()

# instantiate the app
server = Flask(__name__)
server.config.from_object(FlaskConfig())
Bootstrap(server)

app = dash.Dash(server=server, show_undo_redo=False, url_base_pathname="/")
app.title = "MiCV"
app.config.suppress_callback_exceptions = True

# Setup cache
cache = Cache(app.server, config={
    'CACHE_TYPE': 'redis',
    'CACHE_THRESHOLD': 200,  # should be equal to maximum number of active users
    "CACHE_DEFAULT_TIMEOUT": 30000
    }
)

# Setup Flask-Security
user_datastore = SQLAlchemySessionUserDatastore(db_session,
                                                User, Role)
security_ctx = security.init_app(app.server, user_datastore, register_blueprint=True)


def delay_flask_security_mail(msg):
	print("[DEBUG] running delay_flask_security_mail")
	with app.server.app_context():
		send_flask_mail.delay(subject=msg.subject, sender=msg.sender,
		    	              recipients=msg.recipients, body=msg.body,
		    	              html=msg.html)
	return None
security_ctx.send_mail_task(delay_flask_security_mail)

# Create a user to test with
@server.before_first_request
def initialize_database():
	init_db()
	db_session.commit()
'''
def create_test_user():
    init_db()
    try:
    	user_datastore.create_user(email='nigeil@yahoo.com', username="nigeil", password='01010101')
    except exc.SQLAlchemyError:
    	print("[DEBUG] user already exists")
    db_session.commit()
'''