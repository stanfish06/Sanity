#ifndef PYSANITY_SANITY_H
#define PYSANITY_SANITY_H

/**
 * Loop get_gene_expression_level() over all genes with OpenMP.
 * Arrays are flat row-major pointers: a 2D array x[i][j] is stored at x[i * ncol + j].
 * Indices: i = gene (0..num_genes-1), j = cell (0..num_cells-1), k = variance bin (0..num_bins-1).
 *
 * @param counts_gc            [num_genes][num_cells], counts_gc[i][j] -> count of gene i in cell j
 * @param num_genes            number of genes
 * @param num_cells            number of cells
 * @param cell_totals_c        [num_cells], cell_totals_c[j] -> total counts of cell j (i.e. library size)
 * @param gene_totals_g        [num_genes], gene_totals_g[i] -> total counts of gene i; must be > 0,
 *                             so remove zero-count genes before calling
 * @param variance_min         smallest variance bin
 * @param variance_max         largest variance bin
 * @param num_bins             number of log-spaced variance bins; posterior is a mixture of gaussians, one per bin
 * @param variance_method      how v_g is chosen from the bins
 *                               0 MARG: weighted average of the gaussian bins
 *                               1 MLE:  gaussian bin with the highest likelihood
 *                               2 MAP:  highest likelihood / v, an implicit 1/v prior on v_g (default)
 *                               3 EAP:  bin closest to the posterior mean of v_g
 * @param num_threads          OpenMP threads; <= 0 uses the OpenMP default
 * @param verbose              print progress to stderr every 5% of genes
 * @param mean_ltq_g           out [num_genes], mean_ltq_g[i] -> mu of gene i: library-corrected log mean
 *                             expression
 * @param mean_ltq_error_sq_g  out [num_genes], mean_ltq_error_sq_g[i] -> variance of mu of gene i;
 *                             for MLE/MAP/EAP only trigamma(gene_totals_g[i]), no spread across bins
 * @param log_fold_change_gc   out [num_genes][num_cells], log_fold_change_gc[i][j] -> delta of gene i in cell j:
 *                             library-corrected gaussian noise on the log mean expression
 * @param log_fold_change_error_sq_gc  out [num_genes][num_cells], log_fold_change_error_sq_gc[i][j] -> variance of
 *                             delta of gene i in cell j; zero-count cells use an asymmetric error bar
 *                             (SI Eq. 38)
 * @param gene_variance_g      out [num_genes], gene_variance_g[i] -> v_g of gene i: variance of the
 *                             gaussian prior on delta, i.e. the chosen bin
 * @param variance_likelihood_gb  out [num_genes][num_bins], variance_likelihood_gb[i][k] -> posterior probability of
 *                             variance bin k for gene i; may be null to skip
 * @param error_buf            out [error_buf_len], message of the first kernel exception
 * @param error_buf_len        size of error_buf
 * @return 0 on success, 1 if the kernel threw (message in error_buf)
 */
int sanity_fit(const double *counts_gc, int num_genes, int num_cells,
               const double *cell_totals_c, const double *gene_totals_g,
               double variance_min, double variance_max, int num_bins,
               int variance_method, int num_threads, int verbose,
               double *mean_ltq_g, double *mean_ltq_error_sq_g,
               double *log_fold_change_gc, double *log_fold_change_error_sq_gc,
               double *gene_variance_g, double *variance_likelihood_gb,
               char *error_buf, int error_buf_len);

#endif
