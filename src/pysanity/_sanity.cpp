#include "../calc_true_variation_parallel_prior_mu_sigma.cpp"

#include "_sanity.h"

#include <atomic>
#include <cstdio>

// This function is supposed to track main() in calc_true_variation_parallel_prior_mu_sigma.cpp
int sanity_fit(const double *counts_gc, int num_genes, int num_cells,
               const double *cell_totals_ptr, const double *gene_totals_g,
               double variance_min, double variance_max, int num_bins,
               int variance_method, int num_threads, int verbose,
               double *mean_ltq_g, double *mean_ltq_error_sq_g,
               double *log_fold_change_gc, double *log_fold_change_error_sq_gc,
               double *gene_variance_g, double *variance_likelihood_gb,
               char *error_buf, int error_buf_len)
{
    const std::vector<double> cell_totals_c(cell_totals_ptr, cell_totals_ptr + num_cells);
    std::vector<double> log_cell_totals_c(num_cells);
    double sum_cell_totals = 0;
    for (int c = 0; c < num_cells; ++c) {
        log_cell_totals_c[c] = std::log(cell_totals_c[c]);
        sum_cell_totals += cell_totals_c[c];
    }
    double log_variance_step = std::log(variance_max / variance_min) / ((double)num_bins - 1);
    std::vector<double> variance_grid_b(num_bins);
    for (int b = 0; b < num_bins; ++b)
        variance_grid_b[b] = variance_min * std::exp(log_variance_step * b);

    std::atomic<bool> failed(false);
    std::atomic<int> num_done(0);
    std::string failure_message;
    const int report_every = std::max(1, num_genes / 20);

    #pragma omp parallel num_threads(num_threads > 0 ? num_threads : omp_get_max_threads())
    {
        std::vector<double> umi_counts_c(num_cells);
        std::vector<double> fitting_workspace_c(num_cells);
        #pragma omp for schedule(dynamic)
        for (int g = 0; g < num_genes; ++g)
        {
            if (failed.load(std::memory_order_relaxed)) continue;
            try
            {
                const double *row = counts_gc + (size_t)g * num_cells;
                std::copy(row, row + num_cells, umi_counts_c.begin());
                RowComputation result = get_gene_expression_level(umi_counts_c, cell_totals_c, gene_totals_g[g], variance_grid_b, num_cells, num_bins, variance_method, log_cell_totals_c, sum_cell_totals, fitting_workspace_c);
                mean_ltq_g[g] = result.mean_ltq;
                mean_ltq_error_sq_g[g] = result.mean_ltq_error_sq;
                gene_variance_g[g] = result.gene_variance;
                std::copy(result.log_fold_change_c.begin(), result.log_fold_change_c.end(), log_fold_change_gc + (size_t)g * num_cells);
                std::copy(result.log_fold_change_error_sq_c.begin(), result.log_fold_change_error_sq_c.end(), log_fold_change_error_sq_gc + (size_t)g * num_cells);
                if (variance_likelihood_gb)
                    std::copy(result.variance_likelihood_b.begin(), result.variance_likelihood_b.end(), variance_likelihood_gb + (size_t)g * num_bins);
            }
            catch (const std::exception &e)
            {
                #pragma omp critical(pysanity_failure)
                if (!failed.exchange(true)) failure_message = "gene " + std::to_string(g) + ": " + e.what();
            }
            if (verbose)
            {
                int done = ++num_done;
                if (done % report_every == 0 || done == num_genes)
                    std::fprintf(stderr, "Sanity: %d/%d genes processed\n", done, num_genes);
            }
        }
    }

    if (failed)
    {
        std::snprintf(error_buf, error_buf_len, "%s", failure_message.c_str());
        return 1;
    }
    return 0;
}
