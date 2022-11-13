from snakemake import shell

if snakemake.params.ploidy in (1, 2):
    shell(
        "{snakemake.input.consensus} "
        "-c {snakemake.input.cluster} -v {snakemake.input.vcf} "
        "-a {snakemake.input.alt} -r {snakemake.input.ref} "
        "-p {snakemake.params.ploidy} "
        "--soup_out {snakemake.output.amb_rna} --vcf_out {snakemake.output.vcf}"
    )
else:
    with open(snakemake.output.amb_rna, "w") as f:
        f.write("Genotype inference inavailable because the ploidy number " "is not supported.")
    with open(snakemake.output.vcf, "w") as f:
        f.write()
