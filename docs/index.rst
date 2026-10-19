pysanity
========

Python bindings for `Sanity <https://github.com/jmbreda/Sanity>`_
(SAmpling-Noise-corrected Inference of Transcription ActivitY), a Bayesian model that infers
log transcription quotients (LTQs) and their error bars from single-cell UMI counts.

The results should match the ``Sanity`` command-line tool without writing or parsing count tables.

Install
-------

.. code-block:: bash

   pip install sanity-py


Quickstart
----------

With an :class:`~anndata.AnnData` that holds raw counts in a layer:

.. code-block:: python

   import pysanity as ps

   ps.sanity(adata, layer="counts")
   adata.layers["sanity_log_mean"]   # posterior log mean, cells x genes
   adata.layers["sanity_log_var"]    # its variance
   adata.obs["library_size_sanity"]  # total counts per cell

   # posterior predictive samples for a set of cells (Or all cells if cell_idx is None)
   X = ps.sample(adata, cell_idx)

Without AnnData, :func:`pysanity.fit` takes a dense genes x cells count matrix:

.. code-block:: python

   res = ps.fit(counts_gc, variance_method="MAP")
   res.ltq, res.ltq_error_sq          # genes x cells

Every gene and every cell needs at least one count; remove empty ones first.

Variance methods
----------------

Sanity evaluates each gene on a log-spaced grid of variances ``v_g``. ``variance_method`` picks
how the reported values use that grid:

- ``"MAP"`` (default): bin with the highest likelihood / v, an implicit 1/v prior on v_g
- ``"MLE"``: bin with the highest likelihood
- ``"EAP"``: bin closest to the posterior mean of v_g
- ``"MARG"``: weighted average over all bins

Citation
--------

If you use pysanity, please cite the Sanity paper:

Breda J, Zavolan M, van Nimwegen E. Bayesian inference of gene expression states from single-cell
RNA-seq data. *Nature Biotechnology* 39(8), 1008–1016 (2021).
`doi:10.1038/s41587-021-00875-x <https://doi.org/10.1038/s41587-021-00875-x>`_

.. code-block:: bibtex

   @article{breda2021sanity,
     title   = {Bayesian inference of gene expression states from single-cell RNA-seq data},
     author  = {Breda, J{\'e}r{\'e}mie and Zavolan, Mihaela and van Nimwegen, Erik},
     journal = {Nature Biotechnology},
     volume  = {39},
     number  = {8},
     pages   = {1008--1016},
     year    = {2021},
     doi     = {10.1038/s41587-021-00875-x}
   }

.. toctree::
   :maxdepth: 2
   :hidden:

   api
