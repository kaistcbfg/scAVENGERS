# scAVENGERS
[![Documentation Status](https://readthedocs.org/projects/scavengers/badge/?version=latest)](https://scavengers.readthedocs.io/en/latest/?badge=latest)

## Overview
scAVENGERS demultiplexes snATAC-seq data by genotype, referring to the variant information.

## Execution
### 1. Installation
The command below clones the repository and installs dependencies.
```
wget https://github.com/kaistcbfg/scAVENGERS/archive/refs/tags/v1.0.0.tar.gz
tar -xvzf v1.0.0.tar.gz
conda env create -f scAVENGERS/envs/environment.yml
```
### 2. Setting parameters for scAVENGERS pipeline
Via the config file in yaml format, you can set parameters for the execution of the pipeline. The parameters include the path of the input and output data and settings for each program in the pipeline. The format is provided in `config.yaml`. 

### 3. Running scAVENGERS pipeline
Below is a command to run the whole pipeline for demultiplexing. scAVENGERS require an indexed alignment file, an indexed reference genome file, and a line-seperated list of barcode sequences in a text file.
```
conda activate scavengers
$scAVENGERS_directory/scAVENGERS pipeline --configfile config.yaml -j $THREADS
```

### 4. Accessing the result
Running scAVENGERS pipeline results in a tab-delimited file. This result file `clusters.tsv` is structured like the below. To note, the format is compatible with the cluster result file generated from souporcell (https://github.com/wheaton5/souporcell).

|column name|description|
|---|---|
|barcode|barcode sequence|
|status|status of the barcode sequence: singlet or doublet|
|assignment|cluster number where the barcode is assigned. Doublets are expressed in the form {n}/{n}.|
|log_prob_singleton|log singleton probability calculated by troublet|
|log_prob_doublet|log doublet probability calculated by troublet|
|cluster{n}|log likelihood of the assignment on cluster n.|
## Tutorial: demultiplexing human prefrontal cortex scATAC-seq data
First, prepare data by running the script below.
```
wget http://junglab.kaist.ac.kr/Dataset/scAVENGERS_example.tar.gz
tar -xvf scAVENGERS_example.tar.gz
cd scAVENGERS_example
gzip -d genome.fa.gz
```
Then, you can run the whole pipeline by running below:
```
$scAVENGERS_directory/scAVENGERS pipeline -j 10 --configfile config.yaml
```
## Documentation
For further information, please refer to the documentation at https://scavengers.readthedocs.io/en/latest/.
