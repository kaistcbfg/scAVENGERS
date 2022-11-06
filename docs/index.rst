scAVENGERS
==========
scAVENGERS is a pipeline for clustering cell barcodes in multiplexed single cell ATAC-seq sequencing data.
scAVENGERS achieves high performance through two strategies.
- Tolerant variant calling to resolve sparseness of scATAC-seq data
- Probabilistic modelling on how genomic DNA fragments are incorporated into scATAC-seq reads

scAVENGERS modules
---------------------
scAVENGERS provides some modules to aid demultiplexing single cell DNA sequencing data.
- **`scAVENGERS pipeline`**: run variant calling, allele count matrix generation, and demultiplexing at once.
- **`scAVENGERS cluster`**: run demultiplexing script only.

installation and running
---------------------------
1. Clone the repository and install dependencies.
```
git clone https://github.com/kaistcbfg/scAVENGERS
conda env create -f scAVENGERS/envs/environment.yaml
```
2. After you cloned the repository, you can run scAVENGERS by executing the command below.
```
$SCAVENGERS_DIRECTORY/scAVENGERS [pipeline|cluster] [options]
```

citation
-----------
TBD

Future updates
--------------
The feature below may be updated in future.
- Module for assigning clusters to donors, given reference variants
- Supporting heterogenous ploidy for each cell and genomic region
- Developing generative model for doublet and ambient variant detection
- Some additional optimizations

.. toctree::
   :maxdepth: 3
   Running scAVENGERS
      installing scAVENGERS
      excecuting scAVENGERS
      interpreting scAVENGERS results
   scAVENGERS modules
      scAVENGERS pipeline
      scAVENGERS cluster
   Tutorial


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`