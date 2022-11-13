#!/usr/bin/env python

import sys
import argparse
import subprocess
import numpy as np


def get_distance(x, y):
    return np.sum(x == y) / len(x)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("A program to assign donors to demultiplexed clusters")
    parser.add_argument(
        "-g", "--genotype", help="genotype matrix for each donor and locus", required=True,
    )
    parser.add_argument(
        "-i", "--index", help="index of variants selected for clustering", required=True,
    )
    parser.add_argument("-v", "--vcf", help="vcf file used for clustering", required=True)
    parser.add_argument(
        "-r", "--ref_vcf", help="variants for each donor in multisample vcf file", required=True,
    )
    args = parser.parse_args()

    genotypes = np.load(args.genotype)["arr_0"]
    used_variant_indices = np.load(args.index)["arr_0"]
    with subprocess.Popen(["bcftools", "view", args.ref_vcf], stdout=subprocess.PIPE) as proc:
        for row in proc.stdout:
            if row.decode().startswith("##"):
                continue
            donors = row.decode().rstrip("\n").split("\t")[9:]
            if len(donors) != genotypes.shape[0]:
                raise ValueError(
                    f"Number of donors given in genotype matrix {genotypes.shape[0]} does not match the number of donors given in reference vcf file {len(donors)}."
                )
            break

    # Make true variants compatible with genotype matrix
    distances = np.zeros((len(donors), len(genotypes)), dtype="int")
    variant_idx = 0
    with subprocess.Popen(
        ["bedtools", "intersect", "-a", args.vcf, "-b", args.ref_vcf, "-loj", "-iobuf", "4G",],
        stdout=subprocess.PIPE,
        bufsize=1000,
    ) as proc:
        for raw_variant_idx, line in enumerate(proc.stdout):
            if variant_idx >= len(used_variant_indices):
                break
            if raw_variant_idx != used_variant_indices[variant_idx]:
                continue
            contents = line.decode().rstrip("\n").split("\t")
            if contents[10] != ".":
                gts = np.array([content.split(":")[0].count("1") for content in contents[19:]])
                distances += (genotypes[:, variant_idx] == gts.reshape(-1, 1)).T
            variant_idx += 1

    print("# Assignment result")
    for idx, donor in enumerate(donors):
        print(f"## {donor}:")
        print(f"- distances from each cluster: {distances[idx]}")
        print(f"- assigned to cluster {np.argmax(distances[idx])}")
