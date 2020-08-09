from flask_admin.contrib.sqla import ModelView
from flask_security import current_user

class AdminView(ModelView):
	can_create = True
	can_edit = True
	can_delete = True

	def is_accessible(self):
		if ((current_user is None)
		or not (current_user.has_role('admin'))):
			return False
		else:
			return True