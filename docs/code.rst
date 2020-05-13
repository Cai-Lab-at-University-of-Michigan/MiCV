.. _code:

The code
========

MiCV is a python-based web application. It is written in python and calls heavily upon a few key libraries that contributors will want to be somewhat familiar with:

- ``scanpy``, a single-cell RNA-seq analysis package
- ``anndata``, a companion to ``scanpy`` that defines and operates on the main data object representation of a scRNA-seq dataset
- ``dash`` / ``plotly``, a set of libraries that form the core of the web-application interface - these provide all of the plots, dropdown menus, slider bars, tabs, etc.
- ``flask``, a python web server package that goes along with ``dash`` to actually handle web requests


Don't call it a callback
************************

As a GUI based application, much of the code base here is not actually sequencing related but organized around *responding to user inputs* through the use of **callback functions**. Callbacks in this application are functions that have a function decorator attached to them that specifies the following:

- Which user interface elements this function *modifies* (output)
- Which user interface elements this function *responds to directly* (input)
- Which *other* user interface elements need to be polled to respond correctly to the user input (state)

Since each user element item has a unique identifier in this application, we can pass the names of these elements directly into the callback decorators to link elements together with the functions that respond to and update them. Then, each time those elements are interacted with, the callback function will fire and perform the necessary tasks to update their respective elements. For example, the button on the main tab that says "recalculate everything" has a callback function that is "watching" for button presses and waiting to do some calculations and subsequently update the plot below it in response to a button press. 


Code organization
*****************
At the top level, there are 3 files that set up the web application:

- ``app.py``, where the actual dash application is defined and configured, along with the caching provider
- ``layouts.py``, where the page layouts are instantiated
- ``index.py``, which imports from both of the above: this is where the web application is started, the page layouts are actually returned by the web server, and all of the function callbacks are registered to their respective buttons/plots/etc. in the page layouts

For each tab/page in MiCV there is a corresponding folder in the ``src`` directory that contains 4 files:

- ``_components.py``, where buttons, plots, etc. are defined and uniquely identified
- ``_layouts.py``, where the components defined in ``_components.py`` are laid out to fill the tab
- ``_callbacks.py``, where button presses, plot manipulations, etc. from the web front-end are directly responded to and updated by callback functions
- ``_functions.py``, where analytical other processing functions are defined: typically these are called by the callback functions and are where most of the scRNA-seq analysis code actually resides

There is also a ``plotting`` directory that contains all of the code and parameters used to generate specific types of plots, and a ``helper_functions.py`` file that has a mish-mash of functions needed for data loading/saving, caching, and manipulation. These helper functions should be reorganized in the future, and we would love your thoughts on how best to (re)structure that and other aspects of this code base. Most of my day is spent at the bench, so I can't say I'm always keeping up with best practices!