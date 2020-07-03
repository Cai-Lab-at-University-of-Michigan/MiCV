import os
class DefaultConfig(object):
	SECRET_KEY = os.environ.get("SECRET_KEY", b"\xc9j\xa2@k\x04\x0e\x8a\xe9\xb6\xfbA\xdfsU\x05\xdfe\xec@\x05\x0b\xfd\x9a")

	SECURITY_PASSWORD_SALT = '1381971702542417111521178241349271342290814913196'
	SECURITY_REGISTERABLE = True #allows user registration
	DATABASE_URI = 'sqlite:////srv/www/MiCV/databases/misc.db' # user security database