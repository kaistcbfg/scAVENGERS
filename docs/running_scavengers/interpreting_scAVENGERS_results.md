# Interpreting scAVENGERS results
The demultiplexing result after running `scAVENGERS pipeline` or `scAVENGERS cluster` provides the cluster assigned and likelihood to be assigned to each cluster for each barcode.

## Donor assignment
Because scAVENGERS is unsupervised, the clustering result itself is insufficient to discriminate which donor is the barcode from. If you have reference variants for each donor, you can compare the variants generated from scAVENGERS pipeline against these reference variants. The formula below can be a metric to assign donor to cluster.

${Similarity}_{donor}=N({Loci}_{generated} \cap {Loci}_{donor})$

 One thing to be careful is that a signinificant number of variants which can distinguish donors are low-quality variants. Therefore you should not filter out low-quality variants from the variants generated from scAVENGERS pipeline.
 
 This process is implemented in `scAVENGERS assign` program. Running the command below will provide donor assignment for each cluster. Clustering results, a vcf file (`variants.vcf.gz` for example below) generated from multiplexed sample, and a multisample vcf file (`donor_variants.vcf.gz` for example below) containing variants for each donor are required.
 ```
 scAVENGERS cluster -r ref.mtx -a alt.mtx -b barcodes.txt -o $OUTDIR
 scAVENGERS assign -g $OUTDIR/gt_matrix.npz -i $OUTDIR/variant_index.npz -v variants.vcf.gz -r donor_variants.vcf.gz
 ``` 