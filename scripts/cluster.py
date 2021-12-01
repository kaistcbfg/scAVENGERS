import logging
import argparse
import numpy as np
from scipy.io import mmread
from numba import jit
from joblib import Parallel, delayed, cpu_count
from ado_correction import diploid_corrections
from time import time
from tqdm import tqdm

golden_ratio = (5 ** 0.5 + 1) / 2

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
    help="priors for each cluster",
)
parser.add_argument(
    "-d", "--ado_rate", default=0.05, type=float, help="allele dropout rate"
)
parser.add_argument(
    "-s",
    "--stop_criterion",
    default=0.1,
    type=float,
    help="likelihood change to stop the iteration in EM algorithm",
)
parser.add_argument(
    "-t", "--threads", default=1, type=int, help="number of threads"
)
parser.add_argument(
    "-i",
    "--max_iter",
    default=50,
    type=int,
    help="number of maximum iterations for EM algorithm",
)
args = parser.parse_args()


def get_csr_matrix_elements(matrix):
    return matrix.data, matrix.indices, matrix.indptr


def get_counts_from_csr_items(csr_items):
    pass


@jit(nopython=True, cache=True)
def get_binomial_probs(allele_counts, probs):
    comb = 1 + (allele_counts == 1)
    success = probs ** allele_counts
    fail = (1 - probs) ** (2 - allele_counts)

    return comb * success * fail


@jit(nopython=True, cache=True)
def get_corrected_count_prob(ref_count, alt_count, corrections):
    correction_prob, correction_mask = corrections
    return correction_prob[ref_count] * correction_mask[alt_count]


@jit(nopython=True, cache=True)
def get_likelihood(ref_counts, alt_counts, binomial_probs, corrections):
    corrected_count_probs = get_corrected_count_prob(
        ref_counts, alt_counts, corrections
    )
    likelihood = np.prod(
        np.sum((corrected_count_probs.T * binomial_probs), axis=0)
    )
    print(np.log(np.sum((corrected_count_probs.T * binomial_probs), axis=0)))
    # print(likelihood)

    return likelihood


def get_likelihood_matrix(
    ref_barcode_matrix,
    alt_barcode_matrix,
    prob_matrix,
    priors,
    corrections,
    n_jobs,
):
    """
    Summary:
        Get binomial likelihood matrix with allelic drop correction

    Args:
        ref_barcode_matrix:
            tuple of data, indices, indptr of csr matrix for reference 
            read counts of shape (n_barcodes, n_variants)
        alt_barcode_matrix:
            tuple of data, indices, indptr of csr matrix for alternative 
            read counts of shape (n_barcodes, n_variants)
        prob_matrix:
            numpy array of shape (n_clusters, n_variants)
        priors:
            numpy array of shape (n_clusters) containing priors
        corrections:
            function to correct allele dropouts
    
    Returns:
        likelihood_matrix:
            numpy array of shape (n_barcodes, n_clusters)
    """
    n_clusters, n_variants = prob_matrix.shape
    n_barcodes = ref_barcode_matrix.shape[0]
    likelihood_matrix = np.zeros((n_barcodes, n_clusters))

    eff_ref_counts = np.arange(0, 3).reshape(-1, 1)
    binomial_prob_matrix = np.zeros(
        (n_clusters, len(eff_ref_counts), n_variants)
    )
    for idx, probs in enumerate(prob_matrix):
        binomial_prob_matrix[idx] += get_binomial_probs(eff_ref_counts, probs)

    likelihoods = Parallel(n_jobs=n_jobs)(
        delayed(get_likelihood)(
            ref_counts=ref_barcode_matrix[idx1].toarray()[0],
            alt_counts=alt_barcode_matrix[idx1].toarray()[0],
            binomial_probs=binomial_prob_matrix[idx2],
            corrections=corrections,
        )
        for idx1 in tqdm(range(n_barcodes))
        for idx2 in range(n_clusters)
    )

    t = time()
    likelihood_matrix = np.fromiter(
        likelihoods, np.int64, n_barcodes * n_clusters
    ).reshape(n_barcodes, n_clusters)
    print(time() - t)

    return likelihood_matrix * priors


@jit(nopython=True, parallel=True, cache=True)
def get_posterior_matrix(likelihood_matrix, temperature):
    """
    Summary:
        get posteriors in for DAEM algorithm.
    
    Args:
        likelihood_matrix:
            numpy array of shape (n_barcodes, n_clusters) for likelihoods
        temperature:
            Float type temperature parameter for DAEM algorithm.
    
    Returns:
        numpy array of shape (n_barcodes, n_clusters) for posteriors
    """
    annealed_likelihood_matrix = likelihood_matrix ** (1 / temperature)
    posterior_matrix = annealed_likelihood_matrix / np.sum(
        annealed_likelihood_matrix, axis=1
    ).reshape(-1, 1)

    return posterior_matrix.T


@jit(nopython=True, cache=True)
def get_expected_likelihood(corrected_count_probs, posteriors, prob):
    allele_counts = np.array([0, 1, 2])
    binomial_probs = get_binomial_probs(allele_counts, prob)
    log_likelihoods = np.log(binomial_probs @ corrected_count_probs.T)
    weighted_log_likelihoods = log_likelihoods * posteriors
    expected_likelihood = -1 * np.prod(weighted_log_likelihoods)

    return expected_likelihood


@jit(nopython=True, cache=True)
def get_new_prob(corrected_count_probs, posteriors, threshold):
    start, end = 0, 1
    candidate1 = (start + end) / golden_ratio
    candidate2 = end / golden_ratio

    while abs(end - start) >= threshold:
        candidate1_likelihood = get_expected_likelihood(
            corrected_count_probs, posteriors, candidate1
        )
        candidate2_likelihood = get_expected_likelihood(
            corrected_count_probs, posteriors, candidate2
        )
        if candidate1_likelihood < candidate2_likelihood:
            end = candidate2
        else:
            start = candidate1
        candidate1 = (start + end) / golden_ratio
        candidate2 = end / golden_ratio

    return (start + end) / 2


@jit(nopython=True, parallel=True, cache=True)
def get_new_prob_matrix(
    ref_loci_matrix_items,
    alt_loci_matrix_items,
    posterior_matrix,
    threshold,
    corrections,
):
    ref_data, ref_indices, ref_indptr = ref_loci_matrix_items
    alt_data, alt_indices, alt_indptr = alt_loci_matrix_items
    n_clusters, n_barcodes = posterior_matrix.shape
    n_variants = len(ref_indptr) - 1
    prob_matrix = np.zeros((n_clusters, n_variants))

    for idx1 in prange(n_variants):
        ref_counts = np.zeros(n_barcodes, dtype=np.int64)
        alt_counts = np.zeros(n_barcodes, dtype=np.int64)
        ref_begn, ref_end = ref_indptr[idx1], ref_indptr[idx1 + 1]
        alt_begn, alt_end = alt_indptr[idx1], alt_indptr[idx1 + 1]
        ref_counts[ref_indices[ref_begn:ref_end]] += ref_data[ref_begn:ref_end]
        alt_counts[alt_indices[alt_begn:alt_end]] += alt_data[alt_begn:alt_end]
        corrected_count_probs = get_corrected_count_prob(
            ref_counts, alt_counts, corrections
        )

        for idx2 in prange(n_clusters):
            posteriors = posterior_matrix[idx2, :]
            prob_matrix[idx2, idx1] += get_new_prob(
                corrected_count_probs=corrected_count_probs,
                posteriors=posteriors,
                threshold=threshold,
            )

    return prob_matrix


def do_DAEM(
    ref_matrix,
    alt_matrix,
    temperature,
    prob_matrix,
    priors,
    n_iter,
    stop_criterion,
    corrections,
):
    ref_loci_matrix = ref_matrix.tocsr()
    ref_barcode_matrix = ref_loci_matrix.transpose().tocsr()
    alt_loci_matrix = alt_matrix.tocsr()
    alt_barcode_matrix = ref_loci_matrix.transpose().tocsr()

    prev_total_likelihood = None
    for _ in range(n_iter):
        likelihood_matrix = get_likelihood_matrix(
            get_csr_matrix_elements(ref_barcode_matrix),
            get_csr_matrix_elements(alt_barcode_matrix),
            prob_matrix,
            priors,
            corrections,
        )
        total_likelihood = np.prod(np.sum(likelihood_matrix, axis=1))
        posterior_matrix = get_posterior_matrix(likelihood_matrix, temperature)
        prob_matrix = get_new_prob_matrix(
            get_csr_matrix_elements(ref_loci_matrix),
            get_csr_matrix_elements(alt_loci_matrix),
            posterior_matrix,
            stop_criterion,
            corrections,
        )
        total_likelihood = np.prod(np.sum(likelihood_matrix, axis=1))
        if (
            prev_total_likelihood is not None
            and abs(prev_total_likelihood - total_likelihood) <= stop_criterion
        ):
            break
        if temperature <= 1:
            break
        temperature /= 2
        prev_total_likelihood = total_likelihood

    if abs(prev_total_likelihood - total_likelihood) <= stop_criterion:
        logging.warn("The estimates did not converge. Need more iterations.")

    return prob_matrix, likelihood_matrix


# setting number of threads
if args.threads <= cpu_count():
    print(f"{args.threads} cores used")
else:
    raise ValueError(
        f"Not enough cores. {cpu_count()} cores available in maximum."
    )

# import count matrices
logging.info("Importing count matrices...")
ref_matrix = mmread(args.ref).astype(np.int64)
alt_matrix = mmread(args.alt).astype(np.int64)
logging.info("Importing done. %f cells and %f variants" % ref_matrix.shape)

# initialize estimators
logging.info("Calculating likelihoods...")
prob_matrix = np.random.rand(args.clusters, ref_matrix.shape[0])

# initializing temperature
count_matrix = ref_matrix + alt_matrix
init_temperature = np.mean(count_matrix.sum(axis=0))

# get priors
if args.priors is None:
    priors = np.repeat(1 / args.clusters, args.clusters)
elif len(args.priors) == args.clusters:
    priors = np.array(args.priors) / sum(args.priors)
else:
    raise ValueError("Number of priors must equal to the number of clusters.")

# # TEST
# ref_loci_matrix = ref_matrix.tocsr().astype(np.int64)
# ref_barcode_matrix = ref_loci_matrix.transpose().tocsr()
# alt_loci_matrix = alt_matrix.tocsr().astype(np.int64)
# alt_barcode_matrix = ref_loci_matrix.transpose().tocsr()

# posterior_matrix = np.random.rand(args.clusters, ref_barcode_matrix.shape[0])
# prob_matrix = get_new_prob_matrix(
#     get_csr_matrix_elements(ref_loci_matrix),
#     get_csr_matrix_elements(alt_loci_matrix),
#     posterior_matrix,
#     args.stop_criterion,
#     diploid_corrections,
# )

# TEST 2
likelihood_matrix = get_likelihood_matrix(
    ref_matrix.tocsc().transpose(),
    alt_matrix.tocsc().transpose(),
    prob_matrix,
    priors,
    diploid_corrections,
    args.threads,
)
np.savetxt("like.txt", likelihood_matrix, delimiter="\t")
exit()

# get MLEs and likelihoods
mle_prob_matrix, max_likelihood_matrix = do_DAEM(
    ref_matrix,
    alt_matrix,
    init_temperature,
    prob_matrix,
    priors,
    args.max_iter,
    args.stop_criterion,
    diploid_corrections,
)
logging.info("Calculating likelihoods done")

# print likelihoods
np.savetxt("like.txt", max_likelihood_matrix, delimiter="\t")
