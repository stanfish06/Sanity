import numpy as np
from anndata import AnnData
from numba import jit
from scipy.sparse import issparse

from .core import VarianceMethod, fit

NUMERIC_EPSILON: float = 1e-10


def sanity(
    adata: AnnData,
    *,
    layer: str | None = None,
    obsm: str | None = None,
    nbins: int = 160,
    vmin: float = 0.001,
    vmax: float = 50,
    variance_method: VarianceMethod = "MAP",
    n_threads: int | None = None,
    verbose: bool = True,
    key_added_log_mean: str = "sanity_log_mean",
    key_added_log_var: str = "sanity_log_var",
    key_added_library_size: str = "library_size_sanity",
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | dict[str, np.ndarray] | None:
    """Fit the Sanity noise model on raw counts.

    Parameters
    ----------
    adata
        AnnData with raw counts (cells x genes) in ``X``, ``layers[layer]`` or ``obsm[obsm]``.
    layer
        Layer to read counts from. Uses ``adata.X`` if neither ``layer`` nor ``obsm`` is given.
        Do not specify simultaneously with ``obsm``.
    obsm
        ``obsm`` key to read counts from (cells x features). Outputs then go to ``obsm`` instead of ``layers``.
        Do not specify simultaneously with ``layer``.
    nbins
        Number of log-spaced variance bins between ``vmin`` and ``vmax``.
    vmin
        Smallest variance bin.
    vmax
        Largest variance bin.
    variance_method
        How the gene variance v_g is chosen from the bins:

        - ``"MARG"``: weighted average over all bins
        - ``"MLE"``: bin with the highest likelihood
        - ``"MAP"``: bin with the highest likelihood / v, an implicit 1/v prior on v_g
        - ``"EAP"``: bin closest to the posterior mean of v_g
    n_threads
        Number of OpenMP threads. ``None`` uses the OpenMP default.
    verbose
        Print progress to stderr every 5% of genes.
    key_added_log_mean
        Key for the posterior log mean (LTQ) of each cell and gene, in ``layers`` or ``obsm``.
    key_added_log_var
        Key for the posterior variance of the log mean, in ``layers`` or ``obsm``.
    key_added_library_size
        Key in ``obs`` for the total counts of each cell.
    inplace
        Write results into ``adata``. If ``False``, return them as a dict and leave ``adata`` untouched.
    copy
        Write results into a copy of ``adata`` and return it. Cannot be combined with ``inplace=False``.

    Returns
    -------
    ``None`` by default, writing to ``adata``:

    - ``.layers[key_added_log_mean]`` (or ``.obsm``): posterior log mean, cells x genes
    - ``.layers[key_added_log_var]`` (or ``.obsm``): posterior variance of the log mean, cells x genes
    - ``.obs[key_added_library_size]``: total counts per cell

    With ``copy=True``, the modified copy of ``adata``.
    With ``inplace=False``, a dict with keys ``"log_mean"``, ``"log_var"`` and ``"library_size"``.

    Raises
    ------
    ValueError
        If both ``layer`` and ``obsm`` are given, ``copy=True`` is combined with ``inplace=False``,
        the input is empty, or a gene or cell has no counts.
    """
    if layer is not None and obsm is not None:
        raise ValueError("pass either layer or obsm, not both")
    if copy and not inplace:
        raise ValueError("copy=True cannot be used with inplace=False")
    if obsm is not None:
        X = adata.obsm[obsm]
    else:
        X = adata.X if layer is None else adata.layers[layer]
    if X is None:
        raise ValueError("adata.X is empty; pass layer or obsm with the raw counts")
    # TODO: can we make Sanity work with sparse matrix?
    X = X.toarray() if issparse(X) else np.asarray(X)  # ty: ignore[unresolved-attribute]
    # Sanity runs gene-wisely, so transpose to genes x cells
    result = fit(
        X.T,
        variance_method=variance_method,
        vmin=vmin,
        vmax=vmax,
        nbins=nbins,
        n_threads=n_threads,
        verbose=verbose,
    )
    out = {
        "log_mean": result.ltq.T,
        "log_var": result.ltq_error_sq.T,
        "library_size": X.sum(axis=1),
    }
    if not inplace:
        return out

    if copy:
        adata = adata.copy()
    slot = adata.obsm if obsm is not None else adata.layers
    slot[key_added_log_mean] = out["log_mean"]
    slot[key_added_log_var] = out["log_var"]
    adata.obs[key_added_library_size] = out["library_size"]
    return adata if copy else None


@jit(nopython=True, cache=True)
def _sample_posterior_predictive_counts(
    log_mean: np.ndarray,
    log_var: np.ndarray,
    library_size: np.ndarray,
    cell_idx: np.ndarray,
    n_posterior: int = 1000,
    ltq_var_scale: float = 0.1,
) -> np.ndarray:
    n_cells = cell_idx.shape[0]
    n_genes = log_mean.shape[1]
    counts = np.empty((n_cells, n_genes), dtype=np.float64)
    for i in range(n_cells):
        c = cell_idx[i]
        for g in range(n_genes):
            ltq = np.random.normal(
                log_mean[c, g], ltq_var_scale * np.sqrt(log_var[c, g])
            )
            rate = library_size[c] * np.exp(ltq)
            counts[i, g] = np.mean(np.random.poisson(rate, size=n_posterior))
    return counts


def sample(
    adata: AnnData,
    cell_idx: np.ndarray,
    scale_before_pca: bool = False,
    n_pca_comps: int | None = None,
    n_posterior: int = 1000,
    ltq_var_scale: float = 0.1,
    *,
    key_log_mean: str = "sanity_log_mean",
    key_log_var: str = "sanity_log_var",
    key_library_size: str = "library_size_sanity",
    from_obsm: bool = False,
) -> np.ndarray:
    """Sample posterior predictive expression for a set of cells from a fitted Sanity model.

    For each cell and gene, draws a log mean from N(log_mean, (ltq_var_scale * sqrt(log_var))^2),
    averages ``n_posterior`` Poisson draws at rate library_size * exp(log mean), and maps the
    result back to log fractions.

    Parameters
    ----------
    adata
        AnnData after :func:`sanity`, with log mean and log var in ``layers`` (or ``obsm``) and library size in ``obs``.
    cell_idx
        Integer indices of the cells to sample; may repeat, e.g. bootstrap indices.
    scale_before_pca
        Z-score each gene with the mean and std of the fitted log mean across cells.
    n_pca_comps
        If not ``None``, center each gene (unless ``scale_before_pca``) and project onto ``adata.varm["PCs"]``.
        Only acts as a switch: all columns of ``varm["PCs"]`` are used regardless of the value.
    n_posterior
        Number of Poisson draws averaged per cell and gene.
    ltq_var_scale
        Factor on the posterior standard deviation of the log mean.
    key_log_mean
        Key of the posterior log mean in ``layers`` (or ``obsm``).
    key_log_var
        Key of the posterior variance of the log mean in ``layers`` (or ``obsm``).
    key_library_size
        Key of the total counts per cell in ``obs``.
    from_obsm
        Read ``key_log_mean`` and ``key_log_var`` from ``obsm`` instead of ``layers``.

    Returns
    -------
    Array of shape ``len(cell_idx)`` x genes with sampled log fractions,
    or ``len(cell_idx)`` x n_pcs when projected onto ``varm["PCs"]``.

    Raises
    ------
    ValueError
        If ``from_obsm=True`` is combined with ``n_pca_comps``, since ``varm["PCs"]`` is in gene space.
    """
    if from_obsm and n_pca_comps is not None:
        raise ValueError(
            "n_pca_comps projects onto varm['PCs'], which needs gene-space layers, not obsm"
        )
    slot = adata.obsm if from_obsm else adata.layers
    log_mean = np.ascontiguousarray(slot[key_log_mean])
    log_var = np.ascontiguousarray(slot[key_log_var])
    library_size = np.array(adata.obs[key_library_size])
    X = _sample_posterior_predictive_counts(
        log_mean=log_mean,
        log_var=log_var,
        library_size=library_size,
        cell_idx=cell_idx,
        n_posterior=n_posterior,
        ltq_var_scale=ltq_var_scale,
    )
    # map back to Sanity's log-fraction scale
    N_c = library_size[cell_idx]
    X = np.log(X / N_c[:, np.newaxis] + NUMERIC_EPSILON)
    if scale_before_pca:
        X = (X - np.mean(log_mean, axis=0)) / np.std(log_mean, axis=0)
    elif not scale_before_pca and n_pca_comps is not None:
        X = X - np.mean(log_mean, axis=0)
    if n_pca_comps is not None:
        X = X @ adata.varm["PCs"]
    return X
