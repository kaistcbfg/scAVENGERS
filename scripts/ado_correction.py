import numpy as np


# diploid
effective_allele_counts = np.arange(0, 3)
correction_prob = np.array(
    [[1 / 3, 1 / 4, 1 / 3], [1 / 3, 1 / 4, 0], [1 / 3, 0, 0]]
)
correction_mask = np.array([[1, 1, 1], [0, 1, 1], [0, 0, 1]])
diploid_corrections = correction_prob, correction_mask
