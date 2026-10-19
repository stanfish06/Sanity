import numpy as np
from numpy.typing import NDArray

def fit(
    counts: NDArray[np.float64],
    cell_totals: NDArray[np.float64],
    gene_totals: NDArray[np.float64],
    variance_min: float,
    variance_max: float,
    num_bins: int,
    variance_method: int,
    num_threads: int,
    verbose: bool,
    return_likelihood: bool,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64] | None,
]: ...
