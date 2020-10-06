 <img src="https://github.com/Cai-Lab-at-University-of-Michigan/MiCV/blob/master/images/MiCV_logo.png" width="500" height="140">

# MiCV
A Multi-Informatic Cellular Visualization tool

## About
MiCV is a python dash-based web-application that enables researchers to upload raw scRNA-seq data and perform filtering, analysis, and manual annotation. **It is largely an interactive wrapper for functions provided by [scanpy](https://github.com/theislab/scanpy) and [palantir](https://github.com/dpeerlab/Palantir)**.

MiCV is being released in advance of our upcoming publication on type-II neurogenesis in *Drosophila melanogaster*; you can test out this web application at [micv.works](https://micv.works), or read our pre-print at the [biorxiv](https://www.biorxiv.org/content/10.1101/2020.07.02.184549v1).

MiCV is still in a beta state and a few bugs here and there are to be expected. We welcome your bug reports on our issue tracker! 

## Features
MiCV currently supports the following:

### Data types
* Uploading of 10X mapping directory contents
* Uploading of h5ad formatted scanpy data objects

### Filtering
* Down-sampling UMI counts
* Filtering by min/max genes per cell
* Filtering by min cells per gene
* Identifying a user-defined number of highly-variable genes

### Projections/clustering
* Changing k, the number of neighbors in the neighborhood graph (for projections)
* Using [bbknn](https://github.com/Teichlab/bbknn) to do batch corrected neighborhood discovery
* Changing the clustering resolution used by the [louvain/leiden](https://github.com/vtraag/leidenalg) clusting algorithm
* Generating UMAP projections

### Pseudotime
* Calculating a pseudotime trajectory based on a user-selectable starting cell
* Calculating gene expression trends across pseudotime for all genes using PyGAM (no R dependencies)

### Marker gene analysis
* Detection of marker genes for arbitrary combinations of clusters
* Use of any scanpy-implemented test for marker gene significance

### Annotations
* Viewing clustering, pseudotime, and gene expression UMAP projections all at once
* Viewing violin plots of gene expression
* Viewing gene expression trends over pseudotime
* Filtering cells based on any combination of the above plots
* (*Drosophila only*) pulling data from Flybase records for each selected gene

### Saving
* Download your anndata object in h5ad format, ready to load into scanpy
for further analysis outside of the bounds of this app

More features are coming! We welcome your suggestions and pull requests.

## Screenshots
![demo sample](https://github.com/cailabumich/MiCV/blob/master/images/MiCV_sample_demo.gif)

## Installation

### Linux
There are three main components to the MiCV software package:
* A redis caching backend server
* A celery task queue for long-running tasks
* The gunicorn server that handles web requests and "runs" the MiCV server

Redis will need to be installed manually using your distribution's package manager. It does not require any non-default configuration.

We have provided a requirements.txt file that contains a list of all python packages necessary for running MiCV (including celery and gunicorn). Using pip, install these dependencies: `pip install -r requirements.txt`. 

You will need to manually create the following directory on your computer:
`# mkdir -p /srv/www/MiCV/databases`. Then, make sure it is user-writeable using `# chown -R [username_here]:[group_here] /srv/www/MiCV`  

Then, run the redis, celery, and gunicorn servers all together using the provided `./MiCV.sh` script. Point your browser to [the server](http://localhost:8050), and you should be all set to go! 

## Usage tips
* Be patient! Many functions take time (sometimes considerable amounts of it) to process, especially with larger datasets.
* You need not "recalculate everything" when you just want to change the UMAP projection or clustering parameters. Each of those processing sections has a button enabling them to be calculated independently (and thus more rapidly).
* Pseudotime trajectory inference may fail to converge. You might attempt to retry it with the same parameters/starting cell - it often works the second time around!
* Save your anndata object and load it into scanpy directly for access to any and all of scanpy's extended functionality.

## Who to thank?
Many helping hands went into the creation of this web application. It was written and is maintained by [Nigel S. Michki](https://github.com/nigeil) in the [Cai Lab](https://www.cai-lab.org/) at the University of Michigan.

Other intellectual contributors include:
* Logan A. Walker
* Ye Li
* Dawen Cai
* [ your name here ]

This application relies heavily upon the incredible work done by the authors and maintainers of many critical software packages, including:
* [scanpy](https://github.com/theislab/scanpy)
* [palantir](https://github.com/dpeerlab/Palantir)
* [bbknn](https://github.com/Teichlab/bbknn)
* [louvain/leiden](https://github.com/vtraag/leidenalg)
* [dash](https://plot.ly/dash/)
* [docker](https://www.docker.com/)
* [flask](https://flask.palletsprojects.com)
* [python](https://www.python.org/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [seaborn](https://seaborn.pydata.org/)
* And many other dependencies down the line

We thank them for their contributions to open scientific computing and discovery.
