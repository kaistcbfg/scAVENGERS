import argparse
import pandas as pd
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument(
    "-r",
    "--ref",
    required=True,
    type=str,
    help="absolute path of ref allele count matrix",
)
parser.add_argument(
    "-r",
    "--alt",
    required=True,
    type=str,
    help="absolute path of ref allele count matrix",
)
parser.add_argument(
    "-o", "--out", required=True, type=str, help="output directory"
)
args = parser.parse_args()

if not args.out.endswith("/"):
    args.out += "/"

ref_matrix = mmread(args.ref)
alt_matrix = mmread(args.alt)
count_matrix = ref_matrix + alt_matrix

total = pd.Dataframe(count_matrix)
alternative = pd.DataFrame(alt_matrix)

total.to_csv(args.out + "AD.csv")
alternative.to_csv(args.out + "DP.csv")
