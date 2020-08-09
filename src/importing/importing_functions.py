import os
'''

TODO: replace this definition with an import from helper functions
or config module

Currently, best way would be to import from the helper functions module, 
however it leads to a circular import of the cache function which leads
to failure to load the app.

'''
user_dataset_path  = "/srv/www/MiCV/user_datasets/"

def get_dataset_list(username):
    dataset_path = user_dataset_path + str(username) + "/"
    if (os.path.exists(dataset_path)):
        files = os.listdir(dataset_path)
        if (len(files) > 0):
            files = [[f.replace(".zarr", ""), os.path.getmtime(dataset_path + str(f))] for f in files]
            files.sort(key=lambda x: x[1], reverse=True) # sort by time
            return files
    return []
