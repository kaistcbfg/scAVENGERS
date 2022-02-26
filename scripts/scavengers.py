#!/usr/bin/env python

import os
import sys
import subprocess
import argparse


parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--bam",
    required=True,
    type=str,
    help="alignment file of multiplexed scATAC-seq data",
)
parser.add_argument(
    "-f", "--fasta", required=True, type=str, help="reference genome file",
)
parser.add_argument(
    "-b",
    "--barcodes",
    required=True,
    type=str,
    help="line-seperated list of barcode sequences",
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    type=str,
    help="directory where the output files are stored",
)
parser.add_argument(
    "-k", "--clusters", required=True, type=int, help="number of genotypes"
)
parser.add_argument(
    "--priors",
    required=False,
    type=float,
    nargs="+",
    help="number or proportion of cells in each genotype. Defaults to 1/k.",
)
parser.add_argument(
    "--ploidy", default=2, type=int, help="ploidy. Defaults to 2."
)
parser.add_argument(
    "--err_rate",
    default=0.001,
    type=float,
    help="rate of erraneous variant assignment on a cell barcode",
)
parser.add_argument(
    "--doublet_rate",
    default=0.2,
    type=float,
    help="Maximum difference of normalized likelihood to detect doublets",
)
parser.add_argument(
    "--varq", default=100, type=float, help="Minimum variant quality"
)
parser.add_argument(
    "--read_count", default=10, type=float, help="Minimum read count on loci"
)
parser.add_argument(
    "--stop_criterion",
    default=0.1,
    type=float,
    help="log likelihood change to define convergence for EM algorithm",
)
parser.add_argument(
    "--max_iter",
    default=1000,
    type=int,
    help="number of maximum iterations for a temperature step",
)
parser.add_argument(
    "-t", "--threads", default=1, type=int, help="number of threads"
)
parser.add_argument("-v", "--vcf", type=str, help="Vcf file")
parser.add_argument(
    "-r", "--ref", type=str, help="Reference allele count matrix in mtx format"
)
parser.add_argument(
    "-a", "--alt", type=str, help="Alternate allele count matrix in mtx format"
)
args = parser.parse_args()
args.bam = os.path.abspath(args.bam)
args.fasta = os.path.abspath(args.fasta)
if not os.path.exists(args.output):
    os.mkdir(args.output)
args.output = os.path.abspath(args.output)

# call variants with quality filter
if os.path.exists(f"{args.output}/variants.vcf.gz") or args.vcf is not None:
    pass
else:
    try:
        cmd = (
            f"freebayes-parallel "
            f"<(fasta_generate_regions.py {args.fasta}.fai 100000) {args.threads} "
            f"-f {args.fasta} -iXu -C 2 -q 30 -n 3 -E 1 -m 30 "
            f"--min-coverage 20 --pooled-continuous --skip-coverage 100000 {args.bam} | "
            f'vcffilter -f "QUAL > {args.varq}" | '
            f"bgzip > {args.output}/variants.vcf.gz"
        )
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        print(proc.stdout.decode())

    except Exception as e:
        print(e)
        if os.path.exists(f"{args.output}/variants.vcf.gz"):
            print("removing vcf file because it may be corrupted...")
            os.remove(f"{args.output}/variants.vcf.gz")
        sys.exit()

# run vartrix
if (
    os.path.exists(f"{args.output}/ref.mtx")
    and os.path.exists(f"{args.output}/alt.mtx")
    or args.ref is not None
    and args.alt is not None
):
    pass
else:
    args.vcf = (
        f"{args.output}/variants.vcf.gz" if args.vcf is None else args.vcf
    )
    try:
        cmd = (
            f"vartrix -b {args.bam} -v {args.vcf} --fasta {args.fasta} "
            f"-c {args.barcodes} "
            f"--out-matrix {args.output}/alt.mtx "
            f"--ref-matrix {args.output}/ref.mtx "
            f"--mapq 30 --scoring-method coverage "
            f"--threads {args.threads}"
        )
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash",
        )
        print(proc.stdout.decode())
    except Exception as e:
        print(e)
        if os.path.exists(f"{args.output}/ref.mtx"):
            print("removing reference mtx file because it may be corrupted...")
            os.remove(f"{args.output}/ref.mtx")
        if os.path.exists(f"{args.output}/alt.mtx"):
            print(
                "removing alternative mtx file because it may be corrupted..."
            )
            os.remove(f"{args.output}/alt.mtx")
        sys.exit()

# cluster cell barcodes by genotype
if os.path.exists(f"{args.output}/clusters_tmp.tsv"):
    pass
else:
    args.ref = f"{args.output}/ref.mtx" if args.ref is None else args.ref
    args.alt = f"{args.output}/alt.mtx" if args.alt is None else args.alt
    try:
        directory = os.path.dirname(os.path.abspath(__file__))
        prior_cmd = (
            f"--priors {' '.join(args.priors)}"
            if args.priors is not None
            else ""
        )
        with open(f"{args.output}/clusters_tmp.tsv", "w") as f:
            proc = subprocess.run(
                (
                    f"python {directory}/cluster.py "
                    f"-r {args.ref} "
                    f"-a {args.alt} "
                    f"-b {args.barcodes} "
                    f"-o {args.output} "
                    f"-k {args.clusters} "
                    f"{prior_cmd} "
                    f"--ploidy {args.ploidy} "
                    f"--err_rate {args.err_rate} "
                    f"--doublet_rate {args.doublet_rate} "
                    f"--stop_criterion {args.stop_criterion} "
                    f"--max_iter {args.max_iter} "
                    f"-t {args.threads}"
                ),
                shell=True,
                executable="/bin/bash",
                stdout=f,
                stderr=subprocess.PIPE,
            )
            if proc.stdout is not None:
                print(proc.stdout.decode())
    except Exception as e:
        print(e)
        if os.path.exists(f"{args.output}/clusters_tmp.tsv"):
            print("removing result tsv file because it may be corrupted...")
            os.remove(f"{args.output}/clusters_tmp.tsv")

# cluster cell barcodes by genotype
if os.path.exists(f"{args.output}/clusters.tsv"):
    pass
else:
    args.ref = f"{args.output}/ref.mtx" if args.ref is None else args.ref
    args.alt = f"{args.output}/alt.mtx" if args.alt is None else args.alt
    try:
        with open(f"{args.output}/clusters.tsv", "w") as f:
            directory = os.path.dirname(
                os.path.dirname(os.path.abspath(__file__))
            )
            proc = subprocess.run(
                [
                    f"{directory}/troublet/target/release/troublet",
                    "--alts",
                    args.alt,
                    "--refs",
                    args.ref,
                    "--clusters",
                    f"{args.output}/clusters_tmp.tsv",
                ],
                stdout=f,
                stderr=subprocess.PIPE,
            )
            if proc.stdout is not None:
                print(proc.stdout.decode())
    except Exception as e:
        print(e)
        if os.path.exists(f"{args.output}/clusters.tsv"):
            print(
                "removing complete result tsv file because it may be corrupted..."
            )
            os.remove(f"{args.output}/clusters.tsv")

