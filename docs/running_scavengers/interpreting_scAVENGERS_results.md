# Interpreting scAVENGERS results
The demultiplexing result after running `scAVENGERS pipeline` or `scAVENGERS cluster` provides the cluster assigned and likelihood to be assigned to each cluster for each barcode.

## Donor assignment
Because scAVENGERS is unsupervised, the clustering result itself is insufficient to discriminate which donor is the barcode from. If you have reference variants for each donor, you can compare the variants generated from scAVENGERS pipeline against these reference variants. The formula below can be a metric to assign donor to cluster.

${Score}_{reference}=\frac{N({Variants}_{generated} \cap {Variants}_{reference})}{N({Variants}_{generated})}$

 One thing to be careful is that a signinificant number of variants which can distinguish donors are low-quality variants. Therefore you should not filter out low-quality variants from the variants generated from scAVENGERS pipeline.