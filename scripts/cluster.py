#!/usr/bin/env python

import argparse
import ctypes
from pickle import dump
import numpy as np
from scipy.io import mmread
from scipy.special import logsumexp
from sklearn.preprocessing import binarize
from numba import jit, vectorize, prange, set_num_threads, get_num_threads
from numba.extending import get_cython_function_address
import pandas as pd
from tqdm import tqdm

addr = get_cython_function_address("scipy.special.cython_special", "binom")
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
comb = functype(addr)


@jit(nopython=True)
def hypergeom(k, N, K, n):
    if N < 0 or K < 0 or K > N or n < 0 or n > N:
        value = 0
    elif k < max(0, n + K - N) or k > min(n, K):
        value = 0
    else:
        value = comb(K, k) * comb(N - K, n - k) / comb(N, n)

    return value


@jit(nopython=True)
def binom(k, n, p):
    if n < 0 or p < 0 or p > 1:
        value = 0
    elif k > n or k < 0:
        value = 0
    else:
        value = comb(n, k) * p ** k * (1 - p) ** (n - k)

    return value


@jit(nopython=True)
def poisson(x, l):
    if l < 0:
        value = 0
    else:
        value = l ** x * np.exp(l) / np.prod(np.arange(x) + 1)

    return value


@jit(nopython=True)
def rescale(a):

    return (a - np.min(a)) / (np.max(a) - np.min(a))


def get_selection_rates(
    count_barcode_matrix, assignment, n_clusters, ploidy, access_rate
):
    selection_rates = np.empty(n_clusters, count_barcode_matrix.shape[1])
    for idx in range(n_clusters):
        selection_rates[idx, :] = (
            (
                1
                - (
                    count_barcode_matrix[
                        np.where(assignment == n_clusters)[0]
                    ].getnnz(axis=1)
                    / count_barcode_matrix.shape[0]
                )
                - access_rate
            )
            / (1 - access_rate)
        ) ** (1 / ploidy)

    return selection_rates


@vectorize("float64(int64, int64, int64, int64, float64, float64, float64)")
def get_log_likelihood(
    ref, alt, real, ploidy, observed_total_prob, selection_rate, base_prob
):
    if ref + alt == 0:
        likelihood = 1
    else:
        likelihood = 0
        for latent_total in range(1, ploidy + 1):
            for latent_alt in range(latent_total + 1):
                observed_alt_prob = binom(
                    alt, ref + alt, latent_alt / latent_total
                )
                latent_alt_prob = hypergeom(
                    latent_alt, ploidy, real, latent_total
                )
                latent_total_prob = binom(latent_total, ploidy, selection_rate)
                likelihood += (
                    observed_alt_prob
                    * observed_total_prob
                    * latent_alt_prob
                    * latent_total_prob
                )
        likelihood = likelihood * (1 - base_prob) + base_prob

    return np.log(likelihood)


@jit(nopython=True, parallel=True)
def get_log_likelihood_matrix(
    ref_barcode_matrix,
    alt_barcode_matrix,
    real_count_matrix,
    total_probs,
    priors,
    ploidy,
    base_prob,
):
    n_clusters = real_count_matrix.shape[0]
    n_barcodes = len(ref_barcode_matrix[1]) - 1
    ref_matrix_data, nz_indptr, total_nz_indices = ref_barcode_matrix
    alt_matrix_data = alt_barcode_matrix[0]
    log_likelihood_matrix = np.zeros((n_barcodes, n_clusters))
    selection_rate = total_probs[0] ** (1 / ploidy)

    for idx1 in prange(n_barcodes):
        ref_cnts = ref_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
        alt_cnts = alt_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
        nz_indices = total_nz_indices[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
        for idx2 in prange(n_clusters):
            real_cnts = real_count_matrix[idx2][nz_indices]
            for idx in range(len(nz_indices)):
                log_likelihood_matrix[idx1, idx2] += get_log_likelihood(
                    ref_cnts[idx],
                    alt_cnts[idx],
                    real_cnts[idx],
                    ploidy,
                    total_probs[ref_cnts[idx] + alt_cnts[idx]],
                    selection_rate,
                    base_prob,
                )
            log_likelihood_matrix[idx1, idx2] += np.log(priors[idx2])

    return log_likelihood_matrix


def get_posterior_matrix(log_likelihood_matrix, temperature):
    annealed_likelihood_matrix = log_likelihood_matrix / temperature
    sum_for_each_cell = logsumexp(annealed_likelihood_matrix, axis=1)
    sum_for_each_cell = sum_for_each_cell.reshape(-1, 1)
    posterior_matrix = np.exp(annealed_likelihood_matrix - sum_for_each_cell).T

    return posterior_matrix


@jit(nopython=True, parallel=True)
def update_real_count_matrix(
    ref_loci_matrix,
    alt_loci_matrix,
    posterior_matrix,
    total_probs,
    ploidy,
    base_prob,
):
    n_clusters = posterior_matrix.shape[0]
    n_variants = ref_loci_matrix[1].shape[0] - 1
    real_count_matrix = np.empty((n_variants, n_clusters), dtype=np.int64)
    ref_matrix_data, nz_indptr, total_nz_indices = ref_loci_matrix
    alt_matrix_data = alt_loci_matrix[0]
    drop_rate = total_probs[0] ** (1 / ploidy)
    nodes = np.arange(ploidy + 1)

    for idx1 in prange(n_variants):
        ref_cnts = ref_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
        alt_cnts = alt_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
        nz_indices = total_nz_indices[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
        log_likelihood_matrix = np.zeros((ploidy + 1, len(nz_indices)))
        for idx in range(len(nz_indices)):
            log_likelihood_matrix[:, idx] += get_log_likelihood(
                ref_cnts[idx],
                alt_cnts[idx],
                nodes,
                ploidy,
                total_probs[ref_cnts[idx] + alt_cnts[idx]],
                drop_rate,
                base_prob,
            )
        for idx2 in prange(n_clusters):
            posteriors = posterior_matrix[idx2][nz_indices]
            max_expected_value = -np.inf
            for idx in np.random.permutation(nodes):
                expected_value = np.dot(log_likelihood_matrix[idx], posteriors)
                if expected_value >= max_expected_value:
                    max_expected_value = expected_value
                    real_count_matrix[idx1, idx2] = nodes[idx]

    return real_count_matrix.T


def get_max_likelihoods(
    ref_matrix,
    alt_matrix,
    real_count_matrix,
    init_temperature,
    priors,
    max_iter,
    stop_criterion,
    total_probs,
    ploidy,
    base_prob,
):
    ref_loci_matrix = ref_matrix
    ref_loci_matrix = (
        ref_loci_matrix.data,
        ref_loci_matrix.indptr,
        ref_loci_matrix.indices,
    )
    ref_barcode_matrix = ref_matrix.tocsc().transpose()
    ref_barcode_matrix = (
        ref_barcode_matrix.data,
        ref_barcode_matrix.indptr,
        ref_barcode_matrix.indices,
    )
    alt_loci_matrix = alt_matrix
    alt_loci_matrix = (
        alt_loci_matrix.data,
        alt_loci_matrix.indptr,
        alt_loci_matrix.indices,
    )
    alt_barcode_matrix = alt_matrix.tocsc().transpose()
    alt_barcode_matrix = (
        alt_barcode_matrix.data,
        alt_barcode_matrix.indptr,
        alt_barcode_matrix.indices,
    )

    temperatures = init_temperature / 2 ** np.arange(np.log2(init_temperature))
    temperatures = np.append(temperatures, [1])
    for temperature in temperatures:
        prev_total_log_likelihood = -np.inf
        for iter in tqdm(range(max_iter)):
            log_likelihood_matrix = get_log_likelihood_matrix(
                ref_barcode_matrix,
                alt_barcode_matrix,
                real_count_matrix,
                total_probs,
                priors,
                ploidy,
                base_prob,
            )
            posterior_matrix = get_posterior_matrix(
                log_likelihood_matrix, temperature
            )
            total_log_likelihood = np.sum(
                log_likelihood_matrix * posterior_matrix.T
            )
            print(f"Iteration {iter}: T={temperature} L={total_log_likelihood}")
            if (
                abs(total_log_likelihood - prev_total_log_likelihood)
                <= stop_criterion
            ):
                break
            real_count_matrix = update_real_count_matrix(
                ref_loci_matrix,
                alt_loci_matrix,
                posterior_matrix,
                total_probs,
                ploidy,
                base_prob,
            )
            prev_total_log_likelihood = total_log_likelihood

    if abs(total_log_likelihood - prev_total_log_likelihood) > stop_criterion:
        print("WARNING: The expected likelihood did not converge.")

    return real_count_matrix, log_likelihood_matrix


# def detect_doublets(log_likelihood_matrix, threshold=None):
# normalized_matrix = log_likelihood_matrix / np.sum(
#     log_likelihood_matrix, axis=1
# ).reshape(-1, 1)
# rescaled_matrix = 1 - rescale(normalized_matrix)
# sorted_rescaled_matrix = np.sort(rescaled_matrix, axis=1)
# difference = np.diff(sorted_rescaled_matrix[:, [-1, -2]], axis=1).flatten()
# doublet_mask = difference < threshold

# return doublet_mask


def cluster(args):
    # set number of threads
    if args.threads <= get_num_threads():
        set_num_threads(args.threads)
        print(f"{args.threads} cores used")
    else:
        raise ValueError(
            f"Not enough cores. {get_num_threads()} cores available in maximum."
        )

    # import count matrices
    ref_matrix = mmread(args.ref).astype(np.int64).tocsr()
    alt_matrix = mmread(args.alt).astype(np.int64).tocsr()
    variant_occurences = binarize(ref_matrix + alt_matrix).sum(axis=1).A1
    variant_indices = np.where(variant_occurences > 10)[0]
    ref_matrix = ref_matrix[variant_indices]
    alt_matrix = alt_matrix[variant_indices]
    print("Importing done. %d variants and %d cell barcodes" % ref_matrix.shape)

    # get priors
    if args.priors is None:
        args.priors = np.repeat(1 / args.clusters, args.clusters)
    elif len(args.priors) == args.clusters:
        args.priors = np.array(args.priors) / sum(args.priors)
    else:
        raise ValueError(
            "Number of priors must equal to the number of clusters."
        )

    # get total read count frequencies
    count_matrix = ref_matrix + alt_matrix
    counts, freqs = np.unique(count_matrix.data, return_counts=True)
    total_freqs = np.zeros(max(counts) + 1)
    total_freqs[counts] += freqs
    total_freqs[0] += count_matrix.shape[0] * count_matrix.shape[1] - sum(freqs)
    total_probs = total_freqs / sum(total_freqs)

    # initialize estimators
    real_count_matrix = np.random.randint(
        0, args.ploidy + 1, (args.clusters, ref_matrix.shape[0])
    )

    # initialize temperature
    init_temperature = np.mean(count_matrix.sum(axis=0)) * 0.1

    # get MLEs and likelihoods
    real_count_matrix, max_likelihood_matrix = get_max_likelihoods(
        ref_matrix=ref_matrix,
        alt_matrix=alt_matrix,
        real_count_matrix=real_count_matrix,
        init_temperature=init_temperature,
        priors=args.priors,
        max_iter=args.max_iter,
        stop_criterion=args.stop_criterion,
        total_probs=total_probs,
        ploidy=args.ploidy,
        base_prob=args.err_rate,
    )
    print("Calculating likelihoods done")

    # write results in a tsv and pkl file
    if not args.output.endswith("/"):
        args.output += "/"
    barcodes = [barcode.strip() for barcode in open(args.barcodes)]
    assignment = np.argmax(max_likelihood_matrix, axis=1).reshape(-1, 1)
    cluster_info = pd.DataFrame(
        np.concatenate(max_likelihood_matrix, axis=1),
        index=barcodes,
    )
    cluster_info.insert(0, "assignment", assignment)
    cluster_info.to_csv(f"{args.output}clusters_tmp.tsv", header=None, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--ref",
        required=True,
        type=str,
        help="Reference allele count matrix in mtx format",
    )
    parser.add_argument(
        "-a",
        "--alt",
        required=True,
        type=str,
        help="Alternate allele count matrix in mtx format",
    )
    parser.add_argument("-v", "--vcf", type=str, help="Vcf file")
    parser.add_argument(
        "-b",
        "--barcodes",
        required=True,
        type=str,
        help="Line-seperated text file of barcode sequences",
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output directory.",
    )
    parser.add_argument(
        "-k", "--clusters", required=True, type=int, help="Number of genotypes."
    )
    parser.add_argument(
        "--priors",
        required=False,
        type=float,
        nargs="+",
        help="Number or proportion of cells in each genotype.",
    )
    parser.add_argument(
        "--ploidy", default=2, type=int, help="Ploidy. Defaults to 2."
    )
    parser.add_argument(
        "--err_rate",
        default=0.001,
        type=float,
        help=(
            "Baseline probability. DO NOT set this parameter zero, because it leads"
            " to log-zeros. Defaults to 0.001."
        ),
    )
    parser.add_argument(
        "--doublet_rate",
        default=0.1,
        type=float,
        help="Maximum difference of normalized likelihood to detect doublets",
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
    args = parser.parse_args()
    cluster(args)
