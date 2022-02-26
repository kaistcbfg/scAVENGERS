# scAVENGERS

## Overview
scAVENGERS demultiplexes snATAC-seq data by genotype, referring to the variant information.

The pipeline is composed of four stages: calling variants, making count matrices for each cell barcode and loci, clustering the cell barcodes, and detecting doublets. The stages are implemented in a Python files scavengers.py and cluster.py, and a precompiled doublet detection program troublet, implemented in Rust.

The generative model is based on data generating process of snATAC-seq reads: the attachment of the reads to barcode sequences and the following PCR amplification of the reads.

## Requirements
scAVENGERS is intended to run at Linux with bash support. The program is successfully tested under the environment below. These requirements are provided in a conda environment or a docker container.
|name|version|
|---|---|
|**External tools**|
|freebayes|1.3.5|
|vartrix|1.1.22|
|**Python and packages**|
|python|3.8.12|
|numba|0.54.1|
|numpy|1.20.3|
|pandas|1.4.0|
|scikit-learn|1.0.2|
|scipy|1.7.3|


## Parameters
To get the information below on command line, use the option `--help`.
|option|abbreviation|description|
|---|---|---|
||||


## Execution
### 1. Installation
#### Docker container
scAVENGERS and its requirements are provided in a docker container.
```
```

#### Running python source code in Anaconda environment
Alternatively, you may download the source code and the requirements manually. The required softwares and packages are provided as a conda environment. Running below will create an environment named 'scavengers.'
```
git clone https://
cd 
conda env create --file environment.yaml
conda activate scavengers
```
### 2. Preparing data
Alignment file in bam format and line-seperated list of barcodes in text file are required to run scAVENGERS.

You can acquire alignment file and barcodes in ease by using Cellranger ATAC count (https://github.com/10XGenomics/cellranger-atac). The tool aligns fastq files by using bwa mem-like algorithm. The alignments are stored in `outs/possorted_bam.bam`, and the barcodes are stored in `outs/filtered_peak_bc_matrix/barcodes.tsv`.

### 3. Running scAVENGERS
The below is a minimal command to run scAVENGERS. With the arguments, the whole pipeline for demultiplexing will be executed. scAVENGERS require an indexed alignment file, an indexed reference genome file and a line-seperated list of barcode sequences in a text file.
```
# docker
docker exec scripts/scavenge.py \
-i {alignment_file.bam} \
-f {reference_genome.fa} \
-b {barcode_file.txt} \
-o {output_directory} \
-k {number_of_genotypes}

#python
python scripts/scavenge.py \
-i {alignment_file.bam} \
-f {reference_genome.fa} \
-b {barcode_file.txt} \
-o {output_directory} \
-k {number_of_genotypes}
```
scAVENGERS will be executed from the middle of the pipeline depending on the arguments `--vcf`, `--ref`, `--alt`.


### 4. Accessing the result
The result file `clusters.tsv` is structured like below. To note, the format is compatitible with the cluster result file generated from souporcell(https://github.com/wheaton5/souporcell).
barcode	status	assignment	log_prob_singleton	log_prob_doublet	cluster0	cluster1	cluster2	cluster3	cluster4

## Tips for better results
### Number of variants
Empirically, the algorithm showed a fair performance for variant count of around a million or two. If the number of variants is excessively low, it is highly recommended to use a bam file generated from another aligner. You may also use another variant caller and skip the variant calling stage by specifing the vcf file by using the option `--vcf`.

```
usage: scavengers.py [-h] -i BAM -f FASTA -b BARCODES -o OUTPUT -k CLUSTERS
                     [--priors PRIORS [PRIORS ...]] [--ploidy PLOIDY]
                     [--err_rate ERR_RATE] [--doublet_rate DOUBLET_RATE]
                     [--varq VARQ] [--read_count READ_COUNT]
                     [--stop_criterion STOP_CRITERION] [--max_iter MAX_ITER]
                     [-t THREADS] [-v VCF] [-r REF] [-a ALT]

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --bam BAM     alignment file of multiplexed scATAC-seq data
  -f FASTA, --fasta FASTA
                        reference genome file
  -b BARCODES, --barcodes BARCODES
                        line-seperated list of barcode sequences
  -o OUTPUT, --output OUTPUT
                        directory where the output files are stored
  -k CLUSTERS, --clusters CLUSTERS
                        number of genotypes
  --priors PRIORS [PRIORS ...]
                        number or proportion of cells in each genotype.
                        Defaults to 1/k.
  --ploidy PLOIDY       ploidy. Defaults to 2.
  --err_rate ERR_RATE   rate of erraneous variant assignment on a cell barcode
  --doublet_rate DOUBLET_RATE
                        Maximum difference of normalized likelihood to detect
                        doublets
  --varq VARQ           Minimum variant quality
  --read_count READ_COUNT
                        Minimum read count on loci
  --stop_criterion STOP_CRITERION
                        log likelihood change to define convergence for EM
                        algorithm
  --max_iter MAX_ITER   number of maximum iterations for a temperature step
  -t THREADS, --threads THREADS
                        number of threads
  -v VCF, --vcf VCF     Vcf file
  -r REF, --ref REF     Reference allele count matrix in mtx format
  -a ALT, --alt ALT     Alternate allele count matrix in mtx format
```
