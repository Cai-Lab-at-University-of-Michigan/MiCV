from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

import zarr
from numcodecs import Blosc
import scipy.sparse as sp

from .celery import task_queue
from layouts import demo
if (demo == False):
    try:
        from configmodule.production_config import ProductionConfig as FlaskConfig
    except:
        from configmodule.default_config import DefaultConfig as FlaskConfig
else:
    from configmodule.default_config import DefaultConfig as FlaskConfig

@task_queue.task
def send_flask_mail(subject=None, sender=None,
                    recipients=None, body=None,
                    html=None):
    print("[DEBUG] running send_flask_mail")
    # Get mail config
    c = FlaskConfig()

    # Construct the message
    message = Mail(
        from_email=sender,
        to_emails=recipients,
        subject=subject,
        plain_text_content=body,
        html_content=html
    )
    try:
        sg = SendGridAPIClient(c.SENDGRID_API_KEY)
        print("[DEBUG] apikey: " + str(c.SENDGRID_API_KEY))
        response = sg.send(message)
        print(response.status_code)
        print(response.body)
        print(response.headers)
    except Exception as e:
        print(e)    
    return None

## Helper zarr caching
@task_queue.task
def write_dense(zarr_cache_dir, key, dense_name, chunk_factors):
    compressor = Blosc(cname='blosclz', clevel=3, shuffle=Blosc.SHUFFLE)
    store = zarr.open(zarr_cache_dir, mode='a')
    
    if (len(store[key]) == 3):
        # assume csr sparse matrix - parse as such
        array_keys = list(store[key].array_keys())
        X = sp.csr_matrix((store[key + "/" + array_keys[0]], 
                                  store[key + "/" + array_keys[1]], store[key + "/" + array_keys[2]]))
    else:
        # assume dense matrix
        # TODO: checking for other cases of sparse matrices/mixed groups
        X = store[key]
    if ((not (dense_name in store))
    or  (X.shape != store[dense_name].shape)) :
        store.create_dataset(dense_name, shape=X.shape,
                             dtype=X.dtype, fill_value=0, 
                             chunks=(int(X.shape[0]/chunk_factors[0]), int(X.shape[1]/chunk_factors[1])),
                             compressor=compressor, overwrite=True)
    if(sp.issparse(X) is True):
        X = X.tocoo()
        store[dense_name].set_coordinate_selection((X.row, X.col), X.data)
    else:
        store[dense_name] = X
    return None
### end celery queue task function definitions