 # scAVENGERS pipeline: running the whole pipeline for demultiplexing scATAC-seq data
`scAVENGERS pipeline` is a wrapper script that runs variant calling, allele count matrix generation, and demultiplexing at once. The pipeline is implemented in snakemake (https://snakemake.github.io/), which is a workflow management system. The pipeline is executed by running the command below.

## Usage
To run scAVENGERS pipeline, you must provide path of the configuration file and number of jobs to use as `--configfile` and `--jobs` options respectively.
```
scAVENGERS pipeline --configfile config.yaml -j $THREADS
```
The final result is given as a tab-seperated file `clusters.tsv` in the output directory. The structure of the result is in the table below.

|column name|description|
|---|---|
|barcode|barcode sequence|
|status|status of the barcode sequence: singlet or doublet|
|assignment|cluster number where the barcode is assigned. Doublets are expressed in a form {n}/{n}.|
|log_prob_singleton|log singleton probability calculated by troublet|
|log_prob_doublet|log doublet probabilty calculated by troublet|
|cluster{n}|log likelihood of the assignment on cluster n.|

## Configuration file
Through configuration file, you can set parameters for each step. The configuration below is provided in `config.yaml` in the scAVENGERS repository.
```
# number of threads to use
# to make this valid, you must set the number of cores when running snakemake pipeline.
# e.g. -j 10
THREADS: 10

# input and output data path
DATA:
  # path to directory for saving output files
  out_dir:
    outdir
  # path to genome sequence file. SHOULD NOT be compressed.
  genome:
    genome.fa
  # path to alignment file
  alignment:
    alignment.bam
  # path to text file containing line-seperated list of barcodes
  barcode:
    barcodes.txt

# cluster.py settings
CLUSTER:
  # number of genotypes in multiplexed sample
  n_genotypes: 10
  # ploidy of the organism of interest
  ploidy: 2
  # error rate correction parameter of probability model
  err_rate: 0.001
  # difference of likelihood to stop the iteration
  stop_criterion: 0.1
  # maximum number of iteration
  max_iter: 1000

# variant caller settings
VARIANT_CALLER:
  # variant caller to use. Choose "freebayes" or "strelka."
  caller: "strelka"
  # minimum quality score of variant, used only if the caller is "freebayes"
  min_qual: 0
  # whether to use low-GQX variants only, used only if the caller is "strelka"
  lowgqx: False

# vartrix settings
VARTRIX:
  # minimum mapping quality of read to use in vartrix
  mapq: 30

# troublet settings
TROUBLET:
  # doublet prior
  doublet_prior: 0.5
  # doublet posterior threshold
  doublet_threshold: 0.9
  # singlet posterior threshold
  singlet_threshold: 0.9
```