import dash_html_components as html

from app import cache

# Progress is an integer percentage in 
# range [0, 100]
def cache_progress(session_ID, progress=None):
	cache_ID = session_ID + "progress"
	if (progress is None):
		ret = cache.get(cache_ID)
	else:
		cache.set(cache_ID, progress)
		ret = progress
	return ret

# History is a list of strings representing the history 
# of operations to the user. Newer items are appended to
# the end of the list
def cache_history(session_ID, history=None):
	cache_ID = session_ID + "history"
	ret = cache.get(cache_ID)
	if (ret is None):
		ret = []
	
	if not (history is None):
		ret.append(history) # add the new item to the status history list
		cache.set(cache_ID, ret)
	
	return ret

def build_history_table(session_ID):
	history = cache_history(session_ID)
	if (history in [None, "", [], 0]):
		return ""
	ret = []
	for item in history[:-10:-1]:
		ret.append(html.P(children=[str(item)]))
	return ret

# State is a dictionary of items representating the 
# state of the adata object.
def cache_state(session_ID, state=None, key=None, val=None):
	cache_ID = session_ID + "state"
	ret = cache.get(cache_ID)
	
	# entire dictionary of state passed - cache it
	if not (state is None):
		cache.set(cache_ID, state)
		ret = state

	# only updating part of the state dictionary
	if (not (key is None)) and (not (val is None)):
		ret[key] = val
		cache.set(cache_ID, ret)
	return ret

# Build the state table up from a state dictionary
# for display using dash_html components
def build_state_table(session_ID):
	state = cache_state(session_ID)
	if (state in ["", [], None, 0]):
		return []
	ret = []

	# Basic information like #cells, #genes, filename, etc.
	keys = [
		["filename", "# cells/obs", "# genes/var", "# counts"], # basic keys
	]

	headers = [
		["Basic information", ""],
	]
    
	for i in range(0, len(headers)):
    	# Append the header for that part of the table
		ret.append(html.Tr([html.Th(col) for col in headers[i]]))

    	# Append the data for that part of the table
		for key in keys[i]:
			if (key in state.keys()):
				ret.append(html.Tr([html.Td(key), html.Td(state[key])]))
			else:
				ret.append(html.Tr([html.Td(key), html.Td("")]))

	return ret
    