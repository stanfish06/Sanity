/* NOTES:
   - this function gets posterior estimates for a single gene
   - mu/var_mu is library-corrected log mean expression
   - delta/var_delta is library-corrected guassian noise for log mean expression (e.g. gaussian noise for rate)
   - variables with _ml are guassian bin with the highest likelhood
    - posterior is a mixture of guassian distributions
    - the regular version computes weighted average of those guassian distributions
   - N_c is the observed libsize
   - n_c is the observed gene counts for a specific gene
*/
void get_gene_expression_level(double *n_c, double *N_c, double n, double vmin, double vmax, double &mu, double &var_mu, double *delta, double *var_delta, int C, int numbin, double a, double b, double *lik, double &v_ml, double &mu_v_ml, double &var_mu_v_ml, double *delta_v_ml, double *var_delta_v_ml);
