import os

OUTDIR = os.path.abspath(config["DATA"]["out_dir"])
SNAKEDIR = os.path.dirname(workflow.snakefile)
STRELKA_DIR = f"{SNAKEDIR}/strelka-2.9.10.centos6_x86_64"
THREADS = config["THREADS"]

rule all:
    input:
        f"{OUTDIR}/clusters.final.tsv",
        f"{OUTDIR}/gt_matrix.npz",
        f"{OUTDIR}/variant_index.npz",
        f"{OUTDIR}/ambient_rna.txt",
        f"{OUTDIR}/cluster_genotypes.vcf"


rule get_external_programs:
    output:
        f"{STRELKA_DIR}/bin/configureStrelkaGermlineWorkflow.py",
        f"{SNAKEDIR}/troublet",
        f"{SNAKEDIR}/consensus.py",
        f"{SNAKEDIR}/stan_consensus.pickle"
    conda: "envs/consensus.yml"
    shell: "{SNAKEDIR}/scripts/get_programs.sh {SNAKEDIR}"


rule call_variants:
    input:
        bam = config["DATA"]["alignment"],
        fasta = config["DATA"]["genome"],
        strelka = f"{STRELKA_DIR}/bin/configureStrelkaGermlineWorkflow.py"
    output:
        f"{OUTDIR}/variants.vcf.gz"
    params:
        min_qual = config["VARIANT_CALLER"]["min_qual"],
    threads:
        THREADS
    run:
        import os
        
        shell("mkdir -p {OUTDIR}")
        if config["VARIANT_CALLER"]["caller"] == "freebayes":
            shell(
                "freebayes-parallel "
                "<(fasta_generate_regions.py {input.fasta}.fai 100000) {threads} "
                "-f {input.fasta} -iXu -C 2 -q 30 -n 3 -E 1 -m 30 "
                "--min-coverage 20 --pooled-continuous --skip-coverage 100000 {input.bam} | "
                'vcffilter -f "QUAL > {params.min_qual}" | '
                "vcftools view -Ob -o {output}"
            )
        elif config["VARIANT_CALLER"]["caller"] == "strelka":
            lowgqx_cmd = "-f LowGQX" if config["VARIANT_CALLER"]["lowgqx"] else ""
            if not os.path.exists(f"{OUTDIR}/runWorkflow.py"):    
                shell(
                    "{input.strelka} --bam {input.bam} --referenceFasta {input.fasta} "
                    "--runDir {OUTDIR}"
                )
            shell("{OUTDIR}/runWorkflow.py -m local -j {threads}")
            shell(
                "bcftools view -v snps {lowgqx_cmd} -Ob -o {output} {OUTDIR}/results/variants/variants.vcf.gz"
            )
            shell("rm {OUTDIR}/results/variants/variants.vcf.gz")
        else:
            raise ValueError("Variant caller must be freebayes or strelka.")

rule make_matrix:
    input:
        bam = config["DATA"]["alignment"],
        fasta = config["DATA"]["genome"],
        barcode = config["DATA"]["barcode"],
        vcf = f"{OUTDIR}/variants.vcf.gz"
    output:
        ref = f"{OUTDIR}/ref.mtx",
        alt = f"{OUTDIR}/alt.mtx"
    params:
        mapq = config["VARTRIX"]["mapq"]
    threads:
        THREADS
    shell: "vartrix -b {input.bam} -v {input.vcf} --fasta {input.fasta} -c {input.barcode} "
        + "--out-matrix {output.alt} --ref-matrix {output.ref} "
        + "--mapq {params.mapq} --scoring-method coverage --threads {threads}"

rule make_clusters:
    input:
        ref = f"{OUTDIR}/ref.mtx",
        alt = f"{OUTDIR}/alt.mtx",
        vcf = f"{OUTDIR}/variants.vcf.gz",
        barcode = config["DATA"]["barcode"]
    output:
        f"{OUTDIR}/clusters.tsv",
        f"{OUTDIR}/gt_matrix.npz",
        f"{OUTDIR}/variant_index.npz",
    params:
        k = config["CLUSTER"]["n_genotypes"],
        ploidy = config["CLUSTER"]["ploidy"],
        err_rate = config["CLUSTER"]["err_rate"],
        stop_criterion = config["CLUSTER"]["stop_criterion"],
        max_iter = config["CLUSTER"]["max_iter"]
    threads:
        THREADS
    run:
        shell(
            SNAKEDIR + "/scripts/cluster.py -r {input.ref} -a {input.alt} -b {input.barcode} "
            "-o {OUTDIR} -k {params.k} --ploidy {params.ploidy} "
            "--err_rate {params.err_rate} --stop_criterion {params.stop_criterion} "
            "--max_iter {params.max_iter} -t {threads}"
        )

rule call_doublets:
    input:
        ref = f"{OUTDIR}/ref.mtx",
        alt = f"{OUTDIR}/alt.mtx",
        cluster = f"{OUTDIR}/clusters.tsv",
        troublet = f"{SNAKEDIR}/troublet"
    output:
        f"{OUTDIR}/clusters.final.tsv"
    params:
        doublet_prior = config["TROUBLET"]["doublet_prior"],
        doublet_threshold = config["TROUBLET"]["doublet_threshold"],
        singlet_threshold = config["TROUBLET"]["singlet_threshold"],
    run:
        shell(
            "{input.troublet} --refs {input.ref} --alts {input.alt} "
            "--clusters {input.cluster} --doublet_prior {params.doublet_prior} "
            "--doublet_threshold {params.doublet_threshold} "
            "--singlet_threshold {params.singlet_threshold} > {output}"
        )


rule get_genotype:
    input:
        ref = f"{OUTDIR}/ref.mtx",
        alt = f"{OUTDIR}/alt.mtx",
        vcf = f"{OUTDIR}/variants.vcf.gz",
        cluster=f"{OUTDIR}/clusters.final.tsv",
        consensus=f"{SNAKEDIR}/consensus.py",
        stan=f"{SNAKEDIR}/stan_consensus.pickle"
    output:
        amb_rna="{OUTDIR}/ambient_rna.txt",
        vcf="{OUTDIR}/cluster_genotypes.vcf"
    params:
        ploidy = config["CLUSTER"]["ploidy"]
    conda: "envs/consensus.yml"
    script: "scripts/get_genotype.py"
