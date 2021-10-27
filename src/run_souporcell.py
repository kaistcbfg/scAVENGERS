import os
import argparse
import subprocess

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering."
)
parser.add_argument(
    "-b", "--barcodes", required=True, help="barcodes.tsv from cellranger"
)
parser.add_argument("-v", "--variants", required=True, help="vcf file used")
parser.add_argument(
    "-t", "--threads", required=True, type=int, help="max threads to use"
)
parser.add_argument(
    "-o",
    "--out_dir",
    required=True,
    help="name of directory to place souporcell files",
)
parser.add_argument(
    "-k",
    "--clusters",
    required=True,
    help="number cluster, tbd add easy way to run on a range of k",
)
parser.add_argument(
    "-p",
    "--ploidy",
    required=False,
    default="2",
    help="ploidy, must be 1 or 2, default = 2",
)
parser.add_argument(
    "--min_alt",
    required=False,
    default="10",
    help="min alt to use locus, default = 10.",
)
parser.add_argument(
    "--min_ref",
    required=False,
    default="10",
    help="min ref to use locus, default = 10.",
)
parser.add_argument(
    "--max_loci",
    required=False,
    default="2048",
    help="max loci per cell, affects speed, default = 2048.",
)
parser.add_argument(
    "--restarts",
    required=False,
    default=100,
    type=int,
    help="number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima",
)
parser.add_argument(
    "--common_variants",
    required=False,
    default=None,
    help="common variant loci or known variant loci vcf, must be vs same reference fasta",
)
parser.add_argument(
    "--known_genotypes",
    required=False,
    default=None,
    help="known variants per clone in population vcf mode, must be .vcf right now we dont accept gzip or bcf sorry",
)
parser.add_argument(
    "--known_genotypes_sample_names",
    required=False,
    nargs="+",
    default=None,
    help="which samples in population vcf from known genotypes option represent the donors in your sample",
)
parser.add_argument(
    "--skip_remap",
    required=False,
    default=False,
    type=bool,
    help="don't remap with minimap2 (not recommended unless in conjunction with --common_variants",
)
parser.add_argument(
    "--no_umi",
    required=False,
    default="False",
    help="set to True if your bam has no UMI tag, will ignore/override --umi_tag",
)
parser.add_argument(
    "--umi_tag",
    required=False,
    default="UB",
    help="set if your umi tag is not UB",
)
parser.add_argument(
    "--cell_tag",
    required=False,
    default="CB",
    help="set if your cell barcode tag is not CB",
)
parser.add_argument(
    "--ignore",
    required=False,
    default=False,
    type=bool,
    help="set to True to ignore data error assertions",
)
args = parser.parse_args()


def souporcell(args, ref_mtx, alt_mtx, final_vcf):
    print("running souporcell clustering")
    cluster_file = args.out_dir + "/clusters_tmp.tsv"
    with open(cluster_file, "w") as log:
        with open(args.out_dir + "/clusters.err", "w") as err:
            directory = os.path.dirname(
                os.path.dirname(os.path.realpath(__file__))
            )
            cmd = [
                directory
                + "/tools/souporcell/souporcell/target/release/souporcell",
                "-k",
                args.clusters,
                "-a",
                alt_mtx,
                "-r",
                ref_mtx,
                "--restarts",
                str(args.restarts),
                "-b",
                args.barcodes,
                "--min_ref",
                args.min_ref,
                "--min_alt",
                args.min_alt,
                "--threads",
                str(args.threads),
            ]
            if not (args.known_genotypes == None):
                cmd.extend(["--known_genotypes", final_vcf])
                if not (args.known_genotypes_sample_names == None):
                    cmd.extend(
                        ["--known_genotypes_sample_names"]
                        + args.known_genotypes_sample_names
                    )
            print(" ".join(cmd))
            subprocess.check_call(cmd, stdout=log, stderr=err)
    subprocess.check_call(["touch", args.out_dir + "/clustering.done"])
    return cluster_file


def doublets(args, ref_mtx, alt_mtx, cluster_file):
    print("running souporcell doublet detection")
    doublet_file = args.out_dir + "/clusters.tsv"
    with open(doublet_file, "w") as dub:
        with open(args.out_dir + "/doublets.err", "w") as err:
            directory = os.path.dirname(
                os.path.dirname(os.path.realpath(__file__))
            )
            subprocess.check_call(
                [
                    directory
                    + "/tools/souporcell/troublet/target/release/troublet",
                    "--alts",
                    alt_mtx,
                    "--refs",
                    ref_mtx,
                    "--clusters",
                    cluster_file,
                ],
                stdout=dub,
                stderr=err,
            )
    subprocess.check_call(["touch", args.out_dir + "/troublet.done"])
    return doublet_file


def consensus(args, ref_mtx, alt_mtx, doublet_file):
    print("running co inference of ambient RNA and cluster genotypes")
    directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    subprocess.check_call(
        [
            "python",
            directory + "/tools/souporcell/consensus.py",
            "-c",
            doublet_file,
            "-a",
            alt_mtx,
            "-r",
            ref_mtx,
            "-p",
            args.ploidy,
            "--output_dir",
            args.out_dir,
            "--soup_out",
            args.out_dir + "/ambient_rna.txt",
            "--vcf_out",
            args.out_dir + "/cluster_genotypes.vcf",
            "--vcf",
            args.variants,
        ]
    )
    subprocess.check_call(["touch", args.out_dir + "/consensus.done"])


#### MAIN RUN SCRIPT
ref_mtx = args.out_dir + "/ref.mtx"
alt_mtx = args.out_dir + "/alt.mtx"
if not (os.path.exists(args.out_dir + "/clustering.done")):
    souporcell(args, ref_mtx, alt_mtx, args.variants)
cluster_file = args.out_dir + "/clusters_tmp.tsv"
if not (os.path.exists(args.out_dir + "/troublet.done")):
    doublets(args, ref_mtx, alt_mtx, cluster_file)
doublet_file = args.out_dir + "/clusters.tsv"
if not (os.path.exists(args.out_dir + "/consensus.done")):
    consensus(args, ref_mtx, alt_mtx, doublet_file)
print("done")

#### END MAIN RUN SCRIPT

