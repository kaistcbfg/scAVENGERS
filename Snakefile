rule all:
    input:
        config["DATA"]["out_dir"] + "clusters.tsv"


rule get_external_programs:
    output:
        strelka = "strelka-2.9.10.centos6_x86_64",
        troublet = "/troublet/target/release/troublet"
    shell: "scripts/get_programs.sh"


rule call_variants:
    input:
        bam = config["DATA"]["alignment"],
        fasta = config["DATA"]["genome"],
        strelka = "strelka-2.9.10.centos6_x86_64"
    output:
        config["DATA"]["out_dir"] + "variants.vcf.gz"
    params:
        outdir = config["DATA"]["out_dir"],
        min_qual = config["VARIANT_CALLER"]["min_qual"],
    threads:
        config["THREADS"]
    run:
        shell("mkdir -p {params.outdir}")
        if config["VARIANT_CALLER"]["caller"] == "freebayes":
            shell(
                "freebayes-parallel "
                + "<(fasta_generate_regions.py {input.fasta}.fai 100000) {threads} "
                + "-f {input.fasta} -iXu -C 2 -q 30 -n 3 -E 1 -m 30 "
                + "--min-coverage 20 --pooled-continuous --skip-coverage 100000 {input.bam} | "
                + 'vcffilter -f "QUAL > {params.min_qual}" | '
                + "vcftools view -Ob > {output}"
            )
        elif config["VARIANT_CALLER"]["caller"] == "strelka":
            lowgqx_cmd = "-f LowGQX" if config["VARIANT_CALLER"]["lowgqx"] else ""
            shell(
                "{input.strelka}/bin/configureStrelkaGermlineWorkflow.py "
                "--bam {input.bam} --referenceFasta {input.fasta} --runDir {params.outdir}"
            )
            shell("{params.outdir}/runWorkflow.py -m local -j {threads}")
            shell(
                "bcftools view {params.outdir}/results/variants/variants.vcf.gz | "
                + "awk '{if ($0~\"#\" || length($4)==1 && length($5)==1) {print}}' | "
                + f"bcftools view -Ob {lowgqx_cmd}"
                + "> {output}" 
            )
        else:
            raise ValueError("Variant caller must be freebayes or strelka.")

rule make_matrix:
    input:
        bam = config["DATA"]["alignment"],
        fasta = config["DATA"]["genome"],
        barcode = config["DATA"]["barcode"],
        vcf = config["DATA"]["out_dir"] + "variants.vcf.gz"
    output:
        ref = config["DATA"]["out_dir"] + "ref.mtx",
        alt = config["DATA"]["out_dir"] + "alt.mtx"
    params:
        mapq = config["VARTRIX"]["mapq"]
    threads:
        config["THREADS"]
    shell: "vartrix -b {input.bam} -v {input.vcf} --fasta {input.fasta} -c {input.barcode} "
        + "--out-matrix {output.alt} --ref-matrix {output.ref} "
        + "--mapq {params.mapq} --scoring-method coverage --threads {threads}"

rule make_clusters:
    input:
        ref = config["DATA"]["out_dir"] + "ref.mtx",
        alt = config["DATA"]["out_dir"] + "alt.mtx",
        vcf = config["DATA"]["out_dir"] + "variants.vcf.gz",
        barcode = config["DATA"]["barcode"]
    output:
        config["DATA"]["out_dir"] + "clusters_tmp.tsv"
    params:
        outdir = config["DATA"]["out_dir"],
        k = config["CLUSTER"]["n_genotypes"],
        ploidy = config["CLUSTER"]["ploidy"],
        err_rate = config["CLUSTER"]["err_rate"],
        stop_criterion = config["CLUSTER"]["stop_criterion"],
        max_iter = config["CLUSTER"]["max_iter"]
    threads:
        config["THREADS"]
    shell: "scripts/cluster.py -r {input.ref} -a {input.alt} -b {input.barcode} "
        + "-o {params.outdir} -k {params.k} --ploidy {params.ploidy} "
        + "--err_rate {params.err_rate} --stop_criterion {params.stop_criterion} "
        + "--max_iter {params.max_iter} -t {threads}"

rule call_doublets:
    input:
        ref = config["DATA"]["out_dir"] + "ref.mtx",
        alt = config["DATA"]["out_dir"] + "alt.mtx",
        cluster = config["DATA"]["out_dir"] + "clusters_tmp.tsv",
        troublet = "troublet/target/release/troublet"
    output:
        config["DATA"]["out_dir"] + "clusters.tsv"
    params:
    shell: "{input.troublet} --refs {input.ref} --alts {input.alt} --clusters {input.cluster}"
