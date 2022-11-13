# scAVENGERS assign: assigning cluster results to each donor if the reference genotypyes are given
`scAVENGERS assign` assigns cluster results to each donor. The program counts the number of loci where the number of alternative alleles matches between a pair of donor and cluster. Then, for each donor, a cluster with maximum count of matching loci is assigned.
## Usage
Running the script below will result to `$OUTDIR/gt_matrix`, a matrix of predicted alternative allele counts for each donor and locus, and `$OUTDIR/variant_index.npz`, an array of indices of variants selected for demultiplexing.
```
scAVENGERS cluster -r ref.mtx -a alt.mtx -b barcodes.txt -o $OUTDIR
```
Assuming that the vcf file used for clustering is `variants.vcf.gz` and multisample vcf file containing variants of each donor `donor_variants.vcf.gz`, the script below will print out donor assignment result in stdout.
```
scAVENGERS assign -g $OUTDIR/gt_matrix.npz -i $OUTDIR/variant_index.npz -v variants.vcf.gz -r donor_variants.vcf.gz
```
## Result
The result contains raw counts of overlapping loci for each cluster and locus, and cluster assigned for each donor. For example for a multiplexed sample with 5 donors, the result will be like:
```
# Assignment result
## Donor_0:
- distances from each cluster: [61317 93302 59376 59651 61111]
- assigned to cluster 1
## Donor_1:
- distances from each cluster: [61426 60375 59034 96009 61124]
- assigned to cluster 3
## Donor_2:
- distances from each cluster: [60726 60373 94721 59265 60678]
- assigned to cluster 2
## Donor_3:
- distances from each cluster: [89190 59754 58746 59247 60536]
- assigned to cluster 0
## Donor_4:
- distances from each cluster: [61639 60420 59795 59849 94393]
- assigned to cluster 4
```
## Parameters
```
scAVENGERS/scAVENGERS assign --help
usage: A program to assign donors to demultiplexed clusters [-h] -g GENOTYPE -i INDEX -v VCF -r REF_VCF

optional arguments:
  -h, --help            show this help message and exit
  -g GENOTYPE, --genotype GENOTYPE
                        genotype matrix for each donor and locus
  -i INDEX, --index INDEX
                        index of variants selected for clustering
  -v VCF, --vcf VCF     vcf file used for clustering
  -r REF_VCF, --ref_vcf REF_VCF
                        variants for each donor in multisample vcf file
 ```