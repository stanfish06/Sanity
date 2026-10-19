import numpy as np
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

from pysanity import (
    fit,
    sample,
    sanity,
)


@pytest.fixture
def adata(counts) -> AnnData:
    X = counts.to_numpy()
    X = X[X.sum(axis=1) > 0]
    return AnnData(
        X=np.zeros((X.shape[1], X.shape[0])), layers={"counts": csr_matrix(X.T)}
    )


@pytest.fixture
def expected(adata):
    res = fit(adata.layers["counts"].T.toarray())
    return {
        "log_mean": res.ltq.T,
        "log_var": res.ltq_error_sq.T,
        "library_size": np.asarray(adata.layers["counts"].sum(axis=1)).ravel(),
    }


def test_writes_layers_matching_core(adata, expected):
    assert sanity(adata, layer="counts", verbose=False) is None
    np.testing.assert_allclose(adata.layers["sanity_log_mean"], expected["log_mean"])
    np.testing.assert_allclose(adata.layers["sanity_log_var"], expected["log_var"])
    np.testing.assert_allclose(
        adata.obs["library_size_sanity"], expected["library_size"]
    )


def test_not_inplace_returns_matrices_and_leaves_adata(adata, expected):
    out = sanity(adata, layer="counts", inplace=False, verbose=False)
    assert isinstance(out, dict)
    for key in ("log_mean", "log_var", "library_size"):
        np.testing.assert_allclose(out[key], expected[key])
    assert "sanity_log_mean" not in adata.layers
    assert "library_size_sanity" not in adata.obs


def test_copy_returns_new_adata(adata, expected):
    new = sanity(adata, layer="counts", copy=True, verbose=False)
    assert isinstance(new, AnnData) and new is not adata
    np.testing.assert_allclose(
        np.asarray(new.layers["sanity_log_mean"]), expected["log_mean"]
    )
    assert "sanity_log_mean" not in adata.layers


def test_obsm_input_writes_obsm_with_custom_keys(adata, expected):
    adata.obsm["raw"] = adata.layers["counts"].toarray()
    sanity(
        adata,
        obsm="raw",
        key_added_log_mean="mu",
        key_added_log_var="var",
        key_added_library_size="lib",
        verbose=False,
    )
    np.testing.assert_allclose(adata.obsm["mu"], expected["log_mean"])
    np.testing.assert_allclose(adata.obsm["var"], expected["log_var"])
    np.testing.assert_allclose(adata.obs["lib"], expected["library_size"])
    assert "mu" not in adata.layers

    X = sample(
        adata,
        np.arange(5),
        n_posterior=20,
        key_log_mean="mu",
        key_log_var="var",
        key_library_size="lib",
        from_obsm=True,
    )
    assert X.shape == (5, adata.obsm["raw"].shape[1])


def test_rejects_conflicting_arguments(adata):
    with pytest.raises(ValueError, match="layer or obsm"):
        sanity(adata, layer="counts", obsm="raw")
    with pytest.raises(ValueError, match="inplace=False"):
        sanity(adata, layer="counts", copy=True, inplace=False)


def test_sample(adata):
    sanity(adata, layer="counts", verbose=False)
    cell_idx = np.arange(10)
    X = sample(adata, cell_idx, n_posterior=50)
    assert X.shape == (10, adata.n_vars)
    assert np.all(np.isfinite(X))
