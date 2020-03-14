 <img src="https://github.com/cailabumich/MiCV/blob/master/images/MiCV_logo.png" width="500" height="140">

# MiCV
A Multi-Informatic Cellular Visualization tool

## About
MiCV is a python dash-based web-application that enables researchers to upload raw scRNA-seq data and perform filtering, analysis, and manual annotation. **It is largely an interactive wrapper for functions provided by [scanpy](https://github.com/theislab/scanpy) and [palantir](https://github.com/dpeerlab/Palantir)**.

MiCV is being released in advance of our upcoming publication on type-II neurogenesis in *Drosophila melanogaster*; you can test out this web application at [micv.works](https://micv.works)

MiCV is still in a pre-alpha state and likely full of bugs. We welcome your bug reports on our issue tracker! 

## Features
MiCV currently supports the following:

### Data types
* Uploading of 10X mapping directory contents
* Uploading of h5ad formatted scanpy data objects

### Filtering
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
* Calculating gene expression trends across pseudotime for all genes

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

### The easy-ish way - Docker
A docker container is available for MiCV. Running it will spin up a web server at [localhost:8050/MiCV](http://localhost:8050/MiCV) and enable you to upload (locally) your data for analysis. You can find that docker image [here](https://hub.docker.com/r/nigeil/micv).

#### Linux
Install docker using your system's package manager, then start it using your init system (likely systemd - `systemctl start docker`). You can then pull the docker image using `docker pull nigeil/micv:0.2.1` (replace `0.2.1` with the latest version on the [docker repo](https://hub.docker.com/r/nigeil/micv)) and start it using `docker run --name=MiCV -p 8050:8050 -e DASH_DEBUG=true micv:0.2.1`. Finally, you can access the interface from [localhost:8050/MiCV](http://localhost:8050/MiCV). Status/error messages will show up in your terminal and potentially on the web interface as well, so keep a look out as you're trying to debug.  

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
