import numpy as np
import pytest

from pysanity import VARIANCE_METHODS, fit

RTOL = 1e-10
ATOL = 1e-12
TEXT_ATOL = 5e-7


@pytest.mark.parametrize("method", list(VARIANCE_METHODS))
def test_matches_cli(method, counts, cli_output):
    ref = cli_output(method)
    kept = counts.loc[ref["gene_ids"]]
    res = fit(
        kept.to_numpy(), variance_method=method, n_threads=2, return_likelihood=True
    )

    np.testing.assert_allclose(res.mean_ltq, ref["mean_ltq"], rtol=RTOL, atol=ATOL)
    np.testing.assert_allclose(
        res.gene_variance, ref["gene_variance"], rtol=RTOL, atol=ATOL
    )
    np.testing.assert_allclose(
        res.log_fold_change, ref["log_fold_change"], rtol=RTOL, atol=ATOL
    )
    np.testing.assert_allclose(
        np.sqrt(res.log_fold_change_error_sq),
        ref["log_fold_change_error"],
        rtol=RTOL,
        atol=ATOL,
    )
    np.testing.assert_allclose(
        np.sqrt(res.mean_ltq_error_sq), ref["mean_ltq_error"], rtol=0, atol=TEXT_ATOL
    )
    assert res.variance_likelihood is not None
    np.testing.assert_allclose(
        res.variance_likelihood, ref["variance_likelihood"], rtol=0, atol=TEXT_ATOL
    )


def test_thread_count_does_not_change_results(counts):
    X = counts.to_numpy()
    X = X[X.sum(axis=1) > 0]
    one = fit(X, n_threads=1)
    many = fit(X, n_threads=4)
    np.testing.assert_array_equal(one.log_fold_change, many.log_fold_change)
    np.testing.assert_array_equal(
        one.log_fold_change_error_sq, many.log_fold_change_error_sq
    )


def test_rejects_zero_count_genes(counts):
    with pytest.raises(ValueError, match="genes have no counts"):
        fit(counts.to_numpy())


def test_rejects_unknown_method(counts):
    with pytest.raises(ValueError, match="variance_method"):
        fit(counts.to_numpy(), variance_method="max_v")  # ty: ignore[invalid-argument-type]
