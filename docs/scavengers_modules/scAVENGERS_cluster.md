# scAVENGERS cluster: demultiplexing scATAC-seq data
`scAVENGERS cluster` demultiplexes cell barcodes into each donor in unsupervised manner, given reference and alternate allele count matrices generated by vartrix in mtx format
## Usage
scAVENGERS cluster provides cluster assignment results in a tab-seperated format. 

Running this minimal command will produce results under the directory `$OUTDIR`.
```
scAVENGERS cluster -r ref.mtx -a alt.mtx -b barcodes.txt -o $OUTDIR
```
Because scAVENGERS cluster does not perform doublet detection, excluding doublet barcodes before or after running scAVENGERS cluster is required. To note, the output of scAVENGERS is compatible to troublet in souporcell pipeline, so you can use troublet to detect doublets after demultiplexing.
```
# Strategy 1: Filtering out doublet barcodes after running doublet detection tools
cat clusters.tsv | LC_ALL=C grep -F -f $SINGLET_BARCODES > clusters.final.tsv

# Strategy 2: Using troublet as doublet detection tool
$TROUBLET_DIR/troublet -r ref.mtx -a alt.mtx --clusters clusters.tsv > clusters.final.tsv
```
## Results
### Output files
Under the designated output directory, the three files are generated.
|name|description|
|---|---|
|clusters.tsv|cluster result file|
|gt_matrix.npz|donor-variant matrix of alternative allele counts|
|variant_index.npz|indices of variants used for demultiplexing|
### Structure of cluster result file
The cluster result file `clusters.tsv` contains barcodes, assigned cluster, and likelihood for each cluster in each column.
- The first column contains barcode sequences.
- The second column contains cluster assignment results.
- From the third column, log likelihoods are written.
## Parameters
```
scAVENGERS/scAVENGERS cluster --help
usage: cluster.py [-h] -r REF -a ALT [-v VCF] -b BARCODES -o OUTPUT -k CLUSTERS [--priors PRIORS [PRIORS ...]]
                  [--ploidy PLOIDY] [--coverage COVERAGE] [--err_rate ERR_RATE] [--stop_criterion STOP_CRITERION]
                  [--max_iter MAX_ITER] [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Reference allele count matrix in mtx format
  -a ALT, --alt ALT     Alternate allele count matrix in mtx format
  -v VCF, --vcf VCF     Vcf file
  -b BARCODES, --barcodes BARCODES
                        Line-seperated text file of barcode sequences
  -o OUTPUT, --output OUTPUT
                        Output directory.
  -k CLUSTERS, --clusters CLUSTERS
                        Number of donors.
  --priors PRIORS [PRIORS ...]
                        Number or proportion of cells in each genotype.
  --ploidy PLOIDY       Ploidy. Defaults to 2.
  --coverage COVERAGE   Minimum coverage of variant to use for clustering
  --err_rate ERR_RATE   Baseline probability. DO NOT set this parameter zero, because it leads to log-zeros. Defaults to
                        0.001.
  --stop_criterion STOP_CRITERION
                        log likelihood change to define convergence for EM algorithm
  --max_iter MAX_ITER   number of maximum iterations for a temperature step
  -t THREADS, --threads THREADS
                        number of threads
```
