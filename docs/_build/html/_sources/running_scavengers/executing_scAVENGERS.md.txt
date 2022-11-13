# Executing scAVENGERS

## Running the whole pipeline
### Execution
The whole pipeline is executed by running the command below.
```
conda activate scavengers
$scAVENGERS_DIRECTORY/scAVENGERS pipeline --configfile config.yaml -j $THREADS
```
### Preparing data
Alignment file in bam format and line-seperated list of barcodes in text file are required to run scAVENGERS.

If you are using products from 10X Genomics, you can acquire the alignment file and barcodes in ease by using Cellranger ATAC count (https://github.com/10XGenomics/cellranger-atac). The tool aligns fastq files by using bwa mem-like algorithm.
- The alignments are stored in `outs/possorted_bam.bam`
- The barcodes are stored in `outs/filtered_peak_bc_matrix/barcodes.tsv`.

### Configuring
The configuration file is provided in `config.yaml` in scAVENGERS repository.

## Running the demultiplexing module
The demultiplexing module `scAVENGERS cluster` is executed by running the command below.
```
conda activate scavengers
$scAVENGERS_DIRECTORY/scAVENGERS cluster -a alt.mtx -r ref.mtx -b barcodes.txt -o clusters.tsv
```
### Preparing data
`scAVENGERS cluster` requires three inputs: count matrices for alternative and reference alleles in .mtx format, and a line-seperated list of barcodes in text file.

You can acquire allele count matrices by using VarTrix (https://github.com/10XGenomics/vartrix). Below is a minimal command to get reference and alternative allele counts.
```
vartrix \
-b alignment.bam \
-v variants.vcf.gz \
--fasta reference.fa \
-c barcodes.txt \
--scoring-method coverage \
--out-matrix {output.alt} --ref-matrix {output.ref} 
```

### Configuring
You can set options from the command when you run scAVENGERS cluster.