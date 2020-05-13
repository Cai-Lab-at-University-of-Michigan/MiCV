.. MiCV documentation master file, created by
   sphinx-quickstart on Wed Apr 29 11:07:38 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MiCV's documentation!
================================

.. image:: images/MiCV_logo.png
  :width: 1000
  :alt: MiCV logo


MiCV (a Multi-informatic Cellular Visualization tool) is a web application that enables researchers to upload raw scRNA-seq data and perform filtering, analysis, and manual annotation without the need for any programming experience whatsoever. It is largely an interactive wrapper for functions provided by scanpy and palantir, two phenominal scRNA-seq and pseudotime analysis packages written in python.

If you are new to single-cell RNA-seq analysis techniques, check out our overview guide for a very brief look at what a typical analysis workflow looks like. Then, move on to the getting started guide for information on how to go through that workflow in MiCV with data of your choosing, or dive right in with our "executive summary" guide. 

MiCV is still in an alpha state, meaning that it likely has both some bugs that need to be worked out and some missing features. We welcome your bug reports and feature requests on our issue tracker, and encourage you to consider reading through the code and making contributions!

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   overview
   userguide
   troubleshooting
   contributing
   code