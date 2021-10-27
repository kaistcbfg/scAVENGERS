import os
import argparse
from scipy.io import mmread, mmwrite

parser = argparse.ArgumentParser()
parser.add_argument("--ref_in", required=True)
parser.add_argument("--ref_out", required=True)
parser.add_argument("--alt_in", required=True)
parser.add_argument("--alt_out", required=True)
parser.add_argument("--min_p", default=0)
parser.add_argument("--max_p", default=1)
args = parser.parse_args()
indices = []

with os.popen(f"grep -v '%' {args.ref_in} | head -n 1 | cut -f 2") as f:
    n_barcodes = int(f.read().split(maxsplit=2)[1])

with os.popen(
    f"tail {args.ref_in} -n +4 | grep -v ' 0' | cut -f 1 -d ' ' | sort | uniq -c"
) as f:
    for row in f:
        count, idx = row.split(maxsplit=1)
        if int(count) / n_barcodes >= float(args.min_p) and int(
            count
        ) / n_barcodes <= float(args.max_p):
            indices.append(int(idx) - 1)

ref_matrix = mmread(args.ref_in).tocsr()
alt_matrix = mmread(args.alt_in).tocsr()
mmwrite(args.ref_out, ref_matrix[indices, :], field="integer")
mmwrite(args.alt_out, alt_matrix[indices, :], field="integer")
