import logging
import argparse
from time import time
import numpy as np
from scipy.io import mmread
from scipy.stats import hypergeom, poisson
from numba import jit, prange, config, set_num_threads, get_num_threads
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
#     help="absolute path of line-sepearted barcode sequences",
# )
parser.add_argument(
    "-k", "--clusters", required=True, type=int, help="number of clusters"
)
parser.add_argument(
    "-p",
    "--priors",
    required=False,
    type=float,
    nargs="+",
    help="number/proportion of cells in each genotype",
)
parser.add_argument("-n", "--ploidy", default=2, type=int, help="ploidy")
parser.add_argument(
    "-s",
    "--stop_criterion",
    default=0.1,
    type=float,
    help="log likelihood change defining convergence in EM algorithm",
)
parser.add_argument(
    "-i",
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


def get_poisson_pmf(expected_count, max_count):
    counts = np.arange(max_count + 1)
    pmf = poisson.pmf(counts, expected_count)

    return pmf


@jit(nopython=True, cache=True)
def get_likelihoods(
    ref_counts, alt_counts, eff_counts, alt_count_pmf, poisson_pmf
):
    likelihoods = np.ones(len(ref_counts))
    for idx in range(len(ref_counts)):
        ref, alt, eff = ref_counts[idx], alt_counts[idx], eff_counts[idx]
        if ref + alt >= alt_count_pmf.shape[0]:
            likelihoods[idx] = poisson_pmf[ref + alt]
        else:
            likelihoods[idx] = (
                alt_count_pmf[ref, alt, eff] * poisson_pmf[ref + alt]
            )

    return likelihoods


@jit(nopython=True, cache=True)
def get_log_likelihood(
    ref_counts, alt_counts, eff_counts, n_variants, alt_count_pmf, poisson_pmf
):
    likelihoods = get_likelihoods(
        ref_counts, alt_counts, eff_counts, alt_count_pmf, poisson_pmf
    )
    log_likelihood = np.sum(np.log(likelihoods))
    log_likelihood += (n_variants - len(ref_counts)) * np.log(poisson_pmf[0])

    return log_likelihood


@jit(nopython=True, parallel=True, nogil=True, cache=True, fastmath=True)
def get_log_likelihood_matrix(
    ref_barcode_matrix,
    alt_barcode_matrix,
    eff_count_matrix,
    priors,
    alt_count_pmf,
    poisson_pmf,
):
    n_clusters, n_variants = eff_count_matrix.shape
    n_barcodes = len(ref_barcode_matrix[1])
    log_likelihood_matrix = np.zeros((n_barcodes, n_clusters))
    ref_matrix_data, nz_indptr, total_nz_indices = ref_barcode_matrix
    alt_matrix_data = alt_barcode_matrix[0]

    for idx1 in prange(n_barcodes):
        for idx2 in prange(n_clusters):
            ref_counts = ref_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            alt_counts = alt_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            nz_indices = total_nz_indices[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            eff_counts = eff_count_matrix[idx2][nz_indices]
            log_likelihood_matrix[idx1, idx2] += get_log_likelihood(
                ref_counts=ref_counts,
                alt_counts=alt_counts,
                eff_counts=eff_counts,
                n_variants=n_variants,
                alt_count_pmf=alt_count_pmf,
                poisson_pmf=poisson_pmf,
            )
            log_likelihood_matrix[idx1, idx2] += priors[idx2]

    return log_likelihood_matrix


def get_posterior_matrix(log_likelihood_matrix, temperature):
    correction = (
        np.max((log_likelihood_matrix / temperature) // -700, axis=1) * 700
    )
    annealed_likelihood_matrix = np.exp(
        log_likelihood_matrix / temperature + correction.reshape(-1, 1)
    )
    sum_for_each_cell = np.sum(annealed_likelihood_matrix, axis=1)
    sum_for_each_cell = sum_for_each_cell.reshape(-1, 1)
    posterior_matrix = (annealed_likelihood_matrix / sum_for_each_cell).T

    return posterior_matrix


@jit(nopython=True, cache=True)
def get_expected_log_likelihood(
    ref_counts,
    alt_counts,
    eff_count,
    posteriors,
    implicit_posterior_sum,
    alt_count_pmf,
    poisson_pmf,
):
    likelihoods = get_likelihoods(
        ref_counts,
        alt_counts,
        np.repeat(eff_count, len(ref_counts)),
        alt_count_pmf,
        poisson_pmf,
    )
    expected_log_likelihood = np.sum(posteriors * np.log(likelihoods))
    expected_log_likelihood += implicit_posterior_sum * np.log(poisson_pmf[0])

    return expected_log_likelihood


@jit(nopython=True, cache=True)
def update_eff_count(
    ref_counts,
    alt_counts,
    posteriors,
    posterior_sum,
    alt_count_pmf,
    poisson_pmf,
):
    implicit_posterior_sum = posterior_sum - np.sum(posteriors)
    expected_log_likelihoods = np.zeros(alt_count_pmf.shape[0])
    for eff_count in np.arange(alt_count_pmf.shape[0], dtype=np.int8):
        expected_log_likelihoods[eff_count] = get_expected_log_likelihood(
            ref_counts,
            alt_counts,
            eff_count,
            posteriors,
            implicit_posterior_sum,
            alt_count_pmf,
            poisson_pmf,
        )

    return np.argmax(expected_log_likelihoods)


@jit(nopython=True, parallel=True, nogil=True, cache=True, fastmath=True)
def update_eff_count_matrix(
    ref_loci_matrix,
    alt_loci_matrix,
    posterior_matrix,
    posterior_sums,
    alt_count_pmf,
    poisson_pmf,
):
    n_clusters = posterior_matrix.shape[0]
    n_variants = len(ref_loci_matrix[1])
    eff_count_matrix = np.zeros((n_variants, n_clusters), dtype=np.int8)
    ref_matrix_data, nz_indptr, total_nz_indices = ref_loci_matrix
    alt_matrix_data = alt_loci_matrix[0]

    for idx1 in prange(n_variants):
        for idx2 in prange(n_clusters):
            ref_counts = ref_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            alt_counts = alt_matrix_data[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            nz_indices = total_nz_indices[nz_indptr[idx1] : nz_indptr[idx1 + 1]]
            eff_count_matrix[idx1, idx2] = update_eff_count(
                ref_counts=ref_counts,
                alt_counts=alt_counts,
                posteriors=posterior_matrix[idx2][nz_indices],
                posterior_sum=posterior_sums[idx2],
                alt_count_pmf=alt_count_pmf,
                poisson_pmf=poisson_pmf,
            )

    return eff_count_matrix.T


def perform_daem(
    ref_matrix,
    alt_matrix,
    eff_count_matrix,
    temperature,
    priors,
    max_iter,
    stop_criterion,
    alt_count_pmf,
    poisson_pmf,
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
            eff_count_matrix,
            priors,
            alt_count_pmf,
            poisson_pmf,
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
            temperature /= 2
        if temperature < 1:
            break
        posterior_sums = np.sum(posterior_matrix, axis=1)
        eff_count_matrix = update_eff_count_matrix(
            ref_loci_matrix,
            alt_loci_matrix,
            posterior_matrix,
            posterior_sums,
            alt_count_pmf,
            poisson_pmf,
        )
        prev_total_log_likelihood = total_log_likelihood
    elapsed_time = time() - t
    print(f"{elapsed_time} seconds elapsed for {iter + 1} iterations")
    print(f"average {elapsed_time / (iter + 1)} s/it")

    if abs(total_log_likelihood - prev_total_log_likelihood) > stop_criterion:
        print("WARNING: The expected likelihood did not converge.")

    return eff_count_matrix, log_likelihood_matrix


def main():
    # setting number of threads
    config.THREADING_LAYER = "tbb"
    if args.threads <= get_num_threads():
        set_num_threads(args.threads)
        logging.info(f"{args.threads} cores used")
    else:
        raise ValueError(
            f"Not enough cores. {get_num_threads()} cores available in maximum."
        )

    # import count matrices
    logging.info("Importing count matrices...")
    ref_matrix = mmread(args.ref).astype(np.int8)
    alt_matrix = mmread(args.alt).astype(np.int8)
    logging.info("Importing done. %f cells and %f variants" % ref_matrix.shape)

    # initialize estimators
    logging.info("Calculating likelihoods...")
    eff_count_matrix = np.random.randint(
        0, args.ploidy + 1, size=(args.clusters, ref_matrix.shape[0])
    )

    # initializing temperature
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

    # get pmf of the model and calculate the MSE for total read count fitting
    alt_count_pmf = get_alt_count_pmf(args.ploidy, 1e-3)
    poisson_pmf = get_poisson_pmf(count_matrix.mean(), count_matrix.max())
    # freqs = np.zeros(count_matrix.max() + 1)
    # counts_found, freqs_found = np.unique(count_matrix.data, return_counts=True)
    # implicit_zero_count = (
    #     count_matrix.shape[0] * count_matrix.shape[1] - count_matrix.nnz
    # )
    # freqs[counts_found] += freqs_found
    # freqs[0] += implicit_zero_count
    # poisson_freqs = poisson_pmf * count_matrix.shape[0] * count_matrix.shape[1]
    # mse = np.mean((poisson_freqs - freqs) ** 2)
    # print(f"MSE={mse}")

    # get MLEs and likelihoods
    mle_matrix, max_likelihood_matrix = perform_daem(
        ref_matrix=ref_matrix,
        alt_matrix=alt_matrix,
        eff_count_matrix=eff_count_matrix,
        temperature=init_temperature,
        priors=args.priors,
        max_iter=args.max_iter,
        stop_criterion=args.stop_criterion,
        alt_count_pmf=alt_count_pmf,
        poisson_pmf=poisson_pmf,
    )
    logging.info("Calculating likelihoods done")

    # write results in a tsv file
    clusters = [f"cluster{n}" for n in range(args.clusters)]
    cluster_info = pd.DataFrame(max_likelihood_matrix, columns=clusters)
    count_info = pd.DataFrame(mle_matrix, index=clusters)
    # cluster_info["assignment"] = cluster_info.idxmax(axis=1)
    cluster_info.to_csv("like.tsv", index=False, sep="\t")
    count_info.to_csv("mle.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()
