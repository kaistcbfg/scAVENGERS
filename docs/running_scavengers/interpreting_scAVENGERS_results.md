# Interpreting scAVENGERS results
The demultiplexing result after running `scAVENGERS pipeline` or `scAVENGERS cluster` provides the cluster assigned and likelihood to be assigned to each cluster for each barcode.

## Donor assignment
Because scAVENGERS is unsupervised, the clustering result itself is insufficient to discriminate which donor is the barcode from. If you have reference variants for each donor, you can compare the variants generated from scAVENGERS pipeline against these reference variants. The counts of variants with same GT for each donor and genotype can be used as scores for assignment. One thing to be careful is that a significant number of variants which can distinguish donors are low-quality variants. Therefore, from the variants called from the multiplexed sample, low-quality variants must not be excluded.
 
 This process is implemented in `scAVENGERS assign` program. Running the command below will provide donor assignment for each cluster. Clustering results, a vcf file (`variants.vcf.gz` for example below) generated from multiplexed sample, and a multisample vcf file (`donor_variants.vcf.gz` for example below) containing variants for each donor are required.
 ```
 scAVENGERS cluster -r ref.mtx -a alt.mtx -b barcodes.txt -o $OUTDIR
 scAVENGERS assign -g $OUTDIR/gt_matrix.npz -i $OUTDIR/variant_index.npz -v variants.vcf.gz -r donor_variants.vcf.gz
 ``` 
