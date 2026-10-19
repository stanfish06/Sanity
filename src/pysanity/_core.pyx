# cython: language_level=3, boundscheck=False, wraparound=False
import numpy as np


cdef extern from "_sanity.h" nogil:
    int sanity_fit(const double *counts_gc, int num_genes, int num_cells,
                   const double *cell_totals_c, const double *gene_totals_g,
                   double variance_min, double variance_max, int num_bins,
                   int variance_method, int num_threads, int verbose,
                   double *mean_ltq_g, double *mean_ltq_error_sq_g,
                   double *log_fold_change_gc, double *log_fold_change_error_sq_gc,
                   double *gene_variance_g, double *variance_likelihood_gb,
                   char *error_buf, int error_buf_len)


def fit(
    const double[:, ::1] counts,
    const double[::1] cell_totals,
    const double[::1] gene_totals,
    double variance_min,
    double variance_max,
    int num_bins,
    int variance_method,
    int num_threads,
    bint verbose,
    bint return_likelihood,
):
    cdef int G = counts.shape[0]
    cdef int C = counts.shape[1]
    mean_ltq = np.empty(G, dtype=np.float64)
    mean_ltq_error_sq = np.empty(G, dtype=np.float64)
    log_fold_change = np.empty((G, C), dtype=np.float64)
    log_fold_change_error_sq = np.empty((G, C), dtype=np.float64)
    gene_variance = np.empty(G, dtype=np.float64)
    likelihood = np.empty((G, num_bins), dtype=np.float64) if return_likelihood else None

    cdef double[::1] mu_v = mean_ltq
    cdef double[::1] mu_err_v = mean_ltq_error_sq
    cdef double[:, ::1] lfc_v = log_fold_change
    cdef double[:, ::1] lfc_err_v = log_fold_change_error_sq
    cdef double[::1] var_v = gene_variance
    cdef double[:, ::1] lik_v
    cdef double *lik_ptr = NULL
    if return_likelihood:
        lik_v = likelihood
        lik_ptr = &lik_v[0, 0]

    cdef char error_buf[512]
    cdef int status
    with nogil:
        status = sanity_fit(
            &counts[0, 0], G, C, &cell_totals[0], &gene_totals[0],
            variance_min, variance_max, num_bins, variance_method, num_threads, verbose,
            &mu_v[0], &mu_err_v[0], &lfc_v[0, 0], &lfc_err_v[0, 0], &var_v[0], lik_ptr,
            error_buf, sizeof(error_buf),
        )
    if status != 0:
        raise RuntimeError(f"Sanity failed at {error_buf.decode()}")
    return mean_ltq, mean_ltq_error_sq, log_fold_change, log_fold_change_error_sq, gene_variance, likelihood
