from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _core

VarianceMethod = Literal["MARG", "MLE", "MAP", "EAP"]
VARIANCE_METHODS: dict[str, int] = {"MARG": 0, "MLE": 1, "MAP": 2, "EAP": 3}


@dataclass(frozen=True)
class SanityResult:
    """Per-gene Sanity fit; gene-major arrays (genes x cells)

    Attributes
    ----------
    mean_ltq
        Shape (genes,). mu_g, the mean log transcription quotient (LTQ) of each gene.
    mean_ltq_error_sq
        Shape (genes,). Variance of mu_g.
    log_fold_change
        Shape (genes, cells). delta_gc, log fold change of each cell from the gene's mean LTQ.
    log_fold_change_error_sq
        Shape (genes, cells). Variance of delta_gc.
    gene_variance
        Shape (genes,). v_g, the chosen variance of the log fold changes.
    variance_likelihood
        Shape (genes, bins). Posterior over the variance bins, or ``None`` unless requested.
    """

    mean_ltq: NDArray[np.float64]
    mean_ltq_error_sq: NDArray[np.float64]
    log_fold_change: NDArray[np.float64]
    log_fold_change_error_sq: NDArray[np.float64]
    gene_variance: NDArray[np.float64]
    variance_likelihood: NDArray[np.float64] | None

    @property
    def ltq(self) -> NDArray[np.float64]:
        return self.mean_ltq[:, None] + self.log_fold_change

    @property
    def ltq_error_sq(self) -> NDArray[np.float64]:
        return self.mean_ltq_error_sq[:, None] + self.log_fold_change_error_sq


def fit(
    counts: ArrayLike,
    variance_method: VarianceMethod = "MAP",
    vmin: float = 0.001,
    vmax: float = 50.0,
    nbins: int = 160,
    n_threads: int | None = None,
    verbose: bool = False,
    return_likelihood: bool = False,
) -> SanityResult:
    """Fit Sanity on a dense count matrix.

    Parameters
    ----------
    counts
        Raw counts of shape (genes, cells). Every gene and every cell needs at least one count.
    variance_method
        How the gene variance v_g is chosen from the bins: ``"MARG"``, ``"MLE"``, ``"MAP"`` or ``"EAP"``.
        See :func:`pysanity.sanity`.
    vmin
        Smallest variance bin.
    vmax
        Largest variance bin.
    nbins
        Number of log-spaced variance bins between ``vmin`` and ``vmax``.
    n_threads
        Number of OpenMP threads. ``None`` uses the OpenMP default.
    verbose
        Print progress to stderr every 5% of genes.
    return_likelihood
        Also return the posterior over the variance bins.

    Returns
    -------
    SanityResult
        Per-gene fit. ``.ltq`` and ``.ltq_error_sq`` give the per-cell log mean and its variance.

    Raises
    ------
    ValueError
        If ``variance_method`` is unknown, ``counts`` is not 2D, or a gene or cell has no counts.
    RuntimeError
        If the Sanity kernel fails on a gene.
    """
    if variance_method not in VARIANCE_METHODS:
        raise ValueError(
            f"variance_method must be one of {list(VARIANCE_METHODS)}, got {variance_method!r}"
        )
    X = np.ascontiguousarray(counts, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError(f"counts must be 2D (genes x cells), got shape {X.shape}")
    gene_totals = X.sum(axis=1)
    if np.any(gene_totals <= 0):
        raise ValueError(
            f"{int(np.sum(gene_totals <= 0))} genes have no counts; remove them before running Sanity"
        )
    cell_totals = X.sum(axis=0)
    if np.any(cell_totals <= 0):
        raise ValueError(
            f"{int(np.sum(cell_totals <= 0))} cells have no counts; remove them before running Sanity"
        )
    return SanityResult(
        *_core.fit(
            X,
            cell_totals,
            gene_totals,
            vmin,
            vmax,
            nbins,
            VARIANCE_METHODS[variance_method],
            n_threads or 0,
            verbose,
            return_likelihood,
        )
    )
