# Tutorial: demultiplexing prefrontal cortex scATAC-seq data
We prepared prefrontal cortex scATAC-seq (PRJNA729525) synthetic mixture data for users who are new to scAVENGERS. To note, the clustering performance is inevitably low because we subsampled part of chromosome 1 to reduce the size of the data, so that most of the call barcodes will be unassigned.
## Downloading the data
The data will be prepared after running the script below. The directory `examples` is where the data are in. Note that the genome fasta file must be decompressed. If the genome fasta file is not decompressed, the whole pipeline will halt when running vartrix.
```
wget https://github.com/kaistcbfg/scAVENGERS/releases/download/v0.1.0/examples.tar.gz
tar -xvzf examples.tar.gz
cd examples
gzip -d genome.fa.gz
```
## Executing scAVENGERS
After setting your current directory into `scAVENGERS_example`, running `scAVENGERS pipeline` will produce results under the directory named `outdir`.
```
scAVENGERS pipeline -j 10 --configfile config.yaml
```