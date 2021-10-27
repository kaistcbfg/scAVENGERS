rule all:
    input:
        config["DATA"]["RESULT_PATH"] + "clusters.tsv"

# rule remap_alignments:
#     input:
#         bam = config["DATA"]["ALIGNMENTS"],
#         genome = config["DATA"]["GENOME"]
#     output:
#         config["DATA"]["RESULT_PATH"] + "remapped.bam"
#     shell:
#         ""

rule call_variants:
    input:
        bam = config["DATA"]["ALIGNMENTS"],
        genome = config["DATA"]["GENOME"]
    output:
        config["DATA"]["RESULT_PATH"] + "variants.vcf" 
    shell:
        """tools/freebayes-1.3.4-linux-static-AMD64 \
        -b {input.bam} -f {input.genome} -v {output} """
        + config["SETTINGS"]["freebayes"]

rule make_matrix:
    input:
        bam = config["DATA"]["ALIGNMENTS"],
        genome = config["DATA"]["GENOME"],
        barcode = config["DATA"]["BARCODES"],
        vcf = config["DATA"]["RESULT_PATH"] + "variants.vcf"
    output:
        ref = config["DATA"]["RESULT_PATH"] + "ref.mtx",
        alt = config["DATA"]["RESULT_PATH"] + "alt.mtx"
    shell:
        """tools/vartrix_linux \
        -b {input.bam} -v {input.vcf} --fasta {input.genome} \
        -c {input.barcode} \
        --out-matrix {output.alt} --ref-matrix {output.ref} """
        + config["SETTINGS"]["vartrix"]


rule make_clusters:
    input:
        ref = config["DATA"]["RESULT_PATH"] + "ref.mtx",
        alt = config["DATA"]["RESULT_PATH"] + "alt.mtx",
        vcf = config["DATA"]["RESULT_PATH"] + "variants.vcf",
        barcode = config["DATA"]["BARCODES"]
    output:
        config["DATA"]["RESULT_PATH"] + "clusters.tsv"
    params:
        out_dir = config["DATA"]["RESULT_PATH"]
    shell:
        """python src/run_souporcell.py \
        -b {input.barcode} \
        -v {input.vcf} \
        -o {params.out_dir} """
        + config["SETTINGS"]["souporcell"]
