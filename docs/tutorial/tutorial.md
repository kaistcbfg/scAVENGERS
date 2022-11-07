# Tutorial: demultiplexing prefrontal cortex scATAC-seq data
We prepared prefrontal cortex scATAC-seq (PRJNA729525) synthetic mixture data for users who are new to scAVENGERS. We extracted half of chromosme 1 due to the memory issue. Thus, the clustering performance is inevitably low, so that most of the call barcodes will be unassigned.
## Downloading the data
The data will be prepared after running the script below. The directory `scAVENGERS_example` is where the data are in. Note that the genome fasta file must be decompressed. If the genome fasta file is not decompressed, the whole pipeline will halt when running vartrix.
```
wget http://junglab.kaist.ac.kr/Dataset/scAVENGERS_example.tar.gz
tar -xvf scAVENGERS_example.tar.gz
cd scAVENGERS_example
gzip -d genome.fa.gz
```
### Executing scAVENGERS
After setting your current directory into `scAVENGERS_example`, running `scAVENGERS pipeline` will produce results under the directory named `outdir`.
```
scAVENGERS pipeline -j 10 --configfile config.yaml
```