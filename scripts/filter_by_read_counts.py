import subprocess
import argparse
from numpy import asarray
from scipy.io import mmread, mmwrite

parser = argparse.ArgumentParser()
parser.add_argument(
    "-r",
    "--ref",
    required=True,
    type=str,
    help="absolute path of ref allele count matrix",
)
parser.add_argument(
    "-a",
    "--alt",
    required=True,
    type=str,
    help="absolute path of ref allele count matrix",
)
parser.add_argument(
    "-v", "--vcf", required=True, type=str, help="absolute path of vcf file",
)
parser.add_argument(
    "-o", "--out", required=True, type=str, help="output directory",
)
parser.add_argument(
    "-n", "--depth", default=2, type=int, help="maximum depth of coverage",
)
args = parser.parse_args()

if not args.out.endswith("/"):
    args.out += "/"

ref_matrix = mmread(args.ref).tocsr()
alt_matrix = mmread(args.alt).tocsr()
count_matrix = ref_matrix + alt_matrix

count_filter = (
    asarray(count_matrix.max(axis=1).todense()).flatten() <= args.depth
)
new_ref_matrix = ref_matrix[count_filter, :]
new_alt_matrix = alt_matrix[count_filter, :]

mmwrite(args.out + "ref.mtx", new_ref_matrix, field="integer")
mmwrite(args.out + "alt.mtx", new_alt_matrix, field="integer")
new_vcf = open(args.out + "variants.vcf", "w")
with subprocess.Popen(
    ["bcftools", "view", "-H", args.vcf], stdout=subprocess.PIPE
) as proc:
    for idx, row in enumerate(proc.stdout.readlines()):
        if count_filter[idx]:
            new_vcf.write(row)
new_vcf.close()
subprocess.run(["bgzip", args.out + "variants.vcf"])
subprocess.run(["bcftools", "index", args.out + "variants.vcf"])
