scAVENGERS
==========
scAVENGERS is a pipeline for clustering cell barcodes in multiplexed single cell ATAC-seq sequencing data.
scAVENGERS achieves high performance through two strategies.

* Tolerant variant calling to resolve sparseness of scATAC-seq data

* Probabilistic modelling on how genomic DNA fragments are incorporated into scATAC-seq reads

Running scAVENGERS
------------------

1. Clone the repository and install dependencies.

.. code-block::

   wget https://github.com/kaistcbfg/scAVENGERS/archive/refs/tags/v0.1.0.tar.gz
   tar -xvzf v0.1.0.tar.gz
   conda env create -f scAVENGERS/envs/environment.yaml

2. After you cloned the repository, you can run scAVENGERS by executing the command below.

.. code-block::

   $SCAVENGERS_DIRECTORY/scAVENGERS [pipeline|cluster] [options]

.. toctree::
   :caption: Running scAVENGERS
   :maxdepth: 1
   :hidden:

   running_scavengers/installing_scAVENGERS.md
   running_scavengers/executing_scAVENGERS.md
   running_scavengers/interpreting_scAVENGERS_results.md


scAVENGERS modules
---------------------

scAVENGERS provides some modules to aid demultiplexing single cell DNA sequencing data.

* scAVENGERS pipeline: scAVENGERS pipeline is a whole pipeline for demultiplexing single cell ATAC-seq data.

* scAVENGERS cluster: scAVENGERS cluster is a module for clustering cell barcodes in a multiplexed single cell ATAC-seq data.

.. toctree::
   :caption: scAVENGERS modules
   :maxdepth: 1
   :hidden:

   scavengers_modules/scAVENGERS_pipeline.md
   scavengers_modules/scAVENGERS_cluster.md

.. toctree::
   :caption: Tutorial: demultiplexing synthetic human prefrontal cortex scATAC-seq mixture
   :maxdepth: 1
   :hidden:

   tutorial/tutorial.md


citation
-----------

TBD

Future update plans
-------------------

The feature below may be updated in future.

* Module for assigning clusters to donors, given reference variants

* Supporting heterogenous ploidy for each cell and genomic region

* Developing generative model for doublet and ambient variant detection

* Some additional optimizations
