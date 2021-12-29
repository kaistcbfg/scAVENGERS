#!/usr/bin/env/python3

import argparse
from time import time
from pickle import dump
import numpy as np
from scipy.io import mmread
from scipy.stats import hypergeom, poisson
from scipy.special import logsumexp, comb
from numba import jit, prange, set_num_threads, get_num_threads
import pandas as pd

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
    help="absolute path of alt allele count matrix",
)
# parser.add_argument(
#     "-b",
#     "--barcodes",
#     required=True,
#     type=str,
#     help="absolute path of line-seperated barcode sequences",
# )
parser.add_argument(
    "-o",
    "--output",
    required=True,
    type=str,
    help="directory where the clustering result and genotype files are stored",
)
parser.add_argument(
    "-k", "--clusters", required=True, type=int, help="number of genotypes"
)
parser.add_argument(
    "--priors",
    required=False,
    type=float,
    nargs="+",
    help="number or proportion of cells in each genotype",
)
parser.add_argument(
    "--ploidy", default=2, type=int, help="ploidy. Defaults to 2."
)
parser.add_argument(
    "--err_rate",
    default=0.01,
    type=float,
    help="rate of erraneous variant assignment on a cell barcode",
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
    help="number of maximum iterations for EM algorithm",
)
parser.add_argument(
    "-t", "--threads", default=1, type=int, help="number of threads"
)
args = parser.parse_args()


def get_alt_count_pmf(ploidy, error_rate):
    pmf = np.zeros((ploidy + 1, ploidy + 1, ploidy + 1))
    for ref in range(ploidy + 1):
        for alt in range(ploidy + 1):
            for eff in range(ploidy + 1):
                p_hypergeom = hypergeom.pmf(alt, ploidy, eff, ref + alt)
                pmf[ref, alt, eff] = p_hypergeom * (1 - error_rate) + error_rate

    return pmf


def get_zipoisson_pmf(mean, stdev, max):
    counts = np.arange(max + 1)
    if (stdev ** 2 + mean ** 2 - mean) == 0 or mean - 1 == 0:
        pmf = np.zeros(len(counts))
        pmf[0] += 1
    else:
        zero_rate = (stdev ** 2 - mean) / (stdev ** 2 + mean ** 2 - mean)
        expected_count = (stdev ** 2 + mean ** 2) / mean - 1
        pmf = (1 - zero_rate) * poisson.pmf(counts, expected_count)
        pmf[0] += zero_rate

    return pmf


def get_maximization_unit(ploidy, error_rate):
    unit = np.zeros((ploidy + 1, ploidy + 1, ploidy + 1))
    for ref in range(ploidy + 1):
        for alt in range(ploidy + 1):
            for eff in range(ploidy + 1):
                unit[ref, alt, eff] = comb(ploidy - eff, ref) * comb(eff, alt)

    return unit


@jit(nopython=True, parallel=True)
def get_total_count_freq_matrix(count_loci_matrix):
    count_matrix_data, nz_indptr, n_barcodes = count_loci_matrix
    freq_matrix = np.zeros((len(nz_indptr) - 1, np.max(count_matrix_data) + 1))

    for idx in prange(len(nz_indptr) - 1):
        counts = count_matrix_data[nz_indptr[idx] : nz_indptr[idx + 1]]
        for count in counts:
            freq_matrix[idx, count] += 1
        freq_matrix[idx, 0] += n_barcodes - len(counts)
        freq_matrix[idx] /= n_barcodes

    return freq_matrix


def get_total_count_pmf(freq_matrix):
    pmf = np.zeros(freq_matrix.shape)
    for idx, freqs in enumerate(freq_matrix):
        mean = np.average(np.arange(len(freqs)), weights=freqs)
        stdev = np.sqrt(
            np.average((np.arange(len(freqs)) - mean) ** 2, weights=freqs)
        )
        pmf[idx, :] = get_zipoisson_pmf(mean, stdev, len(freqs) - 1)

    return pmf


@jit(nopython=True, parallel=True, nogil=True, cache=True)
def get_log_likelihood_matrix(
    ref_barcode_matrix,
    alt_barcode_matrix,
    real_count_matrix,
    priors,
    alt_count_pmf,
    total_count_pmf,
):
    n_clusters = real_count_matrix.shape[0]
    n_barcodes = len(ref_barcode_matrix[1]) - 1
    ref_matrix_data, nz_indptr, total_nz_indices = ref_barcode_matrix
    alt_matrix_data = alt_barcode_matrix[0]
    total_count_logpmf_zero = np.sum(np.log(total_count_pmf[:, 0]))
    log_likelihood_matrix = np.empty((n_barcodes, n_clusters))

    for idx1 in prange(n_barcodes):
        for idx2 in prange(n_clusters):
            ref_counts = ref_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            alt_counts = alt_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            nz_indices = total_nz_indices[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            real_counts = real_count_matrix[idx2][nz_indices]
            likelihoods = np.ones(len(ref_counts))
            for idx3 in prange(len(likelihoods)):
                ref, alt, real, locus_index = (
                    ref_counts[idx3],
                    alt_counts[idx3],
                    real_counts[idx3],
                    nz_indices[idx3],
                )
                if total_count_pmf[locus_index, ref + alt] == 0:
                    continue
                if ref + alt == 0:
                    continue
                if ref + alt >= len(alt_count_pmf):
                    likelihoods[idx3] *= total_count_pmf[locus_index, ref + alt]
                else:
                    likelihoods[idx3] *= (
                        alt_count_pmf[ref, alt, real]
                        * total_count_pmf[locus_index, ref + alt]
                    )
            log_likelihood_matrix[idx1, idx2] = np.sum(np.log(likelihoods))
            # log_likelihood_matrix[idx1, idx2] += (
            #     total_count_logpmf_zero
            #     - np.sum(np.log(total_count_pmf[nz_indices, 0]))
            # )
            log_likelihood_matrix[idx1, idx2] += np.log(priors[idx2])

    return log_likelihood_matrix


def get_posterior_matrix(log_likelihood_matrix, temperature):
    annealed_likelihood_matrix = log_likelihood_matrix / temperature
    sum_for_each_cell = logsumexp(annealed_likelihood_matrix, axis=1)
    sum_for_each_cell = sum_for_each_cell.reshape(-1, 1)
    posterior_matrix = np.exp(annealed_likelihood_matrix - sum_for_each_cell).T

    return posterior_matrix


@jit(nopython=True, parallel=True, nogil=True, cache=True)
def update_real_count_matrix(
    ref_loci_matrix,
    alt_loci_matrix,
    posterior_matrix,
    alt_count_pmf,
    total_count_pmf,
):
    n_clusters = posterior_matrix.shape[0]
    n_variants = len(ref_loci_matrix[1]) - 1
    real_count_matrix = np.empty((n_variants, n_clusters), dtype=np.int8)
    ref_matrix_data, nz_indptr, total_nz_indices = ref_loci_matrix
    alt_matrix_data = alt_loci_matrix[0]

    for idx1 in prange(n_variants):
        for idx2 in prange(n_clusters):
            ref_counts = ref_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            alt_counts = alt_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            posteriors = posterior_matrix[idx2][
                total_nz_indices[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            ]
            expected_log_likelihoods = np.zeros(alt_count_pmf.shape[0])
            for real in prange(alt_count_pmf.shape[0]):
                likelihoods = np.ones(len(ref_counts))
                for idx3 in prange(len(likelihoods)):
                    ref, alt = ref_counts[idx3], alt_counts[idx3]
                    if ref + alt >= len(alt_count_pmf) or ref + alt == 0:
                        continue
                    likelihoods[idx3] *= alt_count_pmf[ref, alt, real]
                expected_log_likelihoods[real] = np.sum(
                    (posteriors * np.log(likelihoods))[likelihoods != 1]
                )
            real_count_matrix[idx1, idx2] = np.argmax(expected_log_likelihoods)

    return real_count_matrix.T


def perform_daem(
    ref_matrix,
    alt_matrix,
    real_count_matrix,
    temperature,
    priors,
    max_iter,
    stop_criterion,
    alt_count_pmf,
    total_count_pmf,
):
    ref_loci_matrix = ref_matrix.tocsr()
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
    alt_loci_matrix = alt_matrix.tocsr()
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

    prev_total_log_likelihood = 1
    posterior_matrix = None
    t = time()
    for iter in range(max_iter):
        log_likelihood_matrix = get_log_likelihood_matrix(
            ref_barcode_matrix,
            alt_barcode_matrix,
            real_count_matrix,
            priors,
            alt_count_pmf,
            total_count_pmf,
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
            temperature = int(temperature / 2)
        if temperature < 1:
            break
        real_count_matrix = update_real_count_matrix(
            ref_loci_matrix,
            alt_loci_matrix,
            posterior_matrix,
            alt_count_pmf,
            total_count_pmf,
        )
        prev_total_log_likelihood = total_log_likelihood
    elapsed_time = time() - t
    print(f"{elapsed_time} seconds elapsed for {iter + 1} iterations")
    print(f"average {elapsed_time / (iter + 1)} s/it")

    if abs(total_log_likelihood - prev_total_log_likelihood) > stop_criterion:
        print("WARNING: The expected likelihood did not converge.")

    return real_count_matrix, log_likelihood_matrix


def main():
    # set number of threads
    if args.threads <= get_num_threads():
        set_num_threads(args.threads)
        print(f"{args.threads} cores used")
    else:
        raise ValueError(
            f"Not enough cores. {get_num_threads()} cores available in maximum."
        )

    # import count matrices
    ref_matrix = mmread(args.ref).astype(np.int8)
    alt_matrix = mmread(args.alt).astype(np.int8)
    print("Importing done. %d variants and %d cell barcodes" % ref_matrix.shape)

    # initialize estimators
    real_count_matrix = np.random.randint(
        0, args.ploidy + 1, size=(args.clusters, ref_matrix.shape[0])
    )

    # initialize temperature
    count_matrix = ref_matrix + alt_matrix
    init_temperature = np.mean(count_matrix.sum(axis=0))

    # get priors
    if args.priors is None:
        args.priors = np.repeat(1 / args.clusters, args.clusters)
    elif len(args.priors) == args.clusters:
        args.priors = np.array(args.priors) / sum(args.priors)
    else:
        raise ValueError(
            "Number of priors must equal to the number of clusters."
        )

    # get pmfs for the model
    count_loci_matrix = count_matrix.tocsr()
    alt_count_pmf = get_alt_count_pmf(args.ploidy, args.err_rate)
    total_count_pmf = get_total_count_freq_matrix(
        (
            count_loci_matrix.data,
            count_loci_matrix.indptr,
            count_loci_matrix.shape[1],
        )
    )
    # total_count_pmf = get_total_count_pmf(total_count_freq_matrix)

    # get MLEs and likelihoods
    real_count_matrix, max_likelihood_matrix = perform_daem(
        ref_matrix=ref_matrix,
        alt_matrix=alt_matrix,
        real_count_matrix=real_count_matrix,
        temperature=init_temperature,
        priors=args.priors,
        max_iter=args.max_iter,
        stop_criterion=args.stop_criterion,
        alt_count_pmf=alt_count_pmf,
        total_count_pmf=total_count_pmf,
    )
    print("Calculating likelihoods done")

    # write results in a tsv and pkl file
    if not args.output.endswith("/"):
        args.output += "/"
    clusters = [f"cluster{n}" for n in range(args.clusters)]
    cluster_info = pd.DataFrame(max_likelihood_matrix, columns=clusters)
    assignment = cluster_info.idxmax(axis=1)
    cluster_info["assignment"] = assignment.str.replace("cluster", "")
    cluster_info.to_csv(f"{args.output}results.tsv", index=False, sep="\t")
    with open(f"{args.output}genotypes.pkl", "wb") as f:
        dump(real_count_matrix, f)


if __name__ == "__main__":
    main()
