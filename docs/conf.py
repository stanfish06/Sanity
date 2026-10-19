from importlib.metadata import version as _version

project = "pysanity"
author = "Zhiyuan(Stan) Yu"
release = _version("sanity-py")
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

napoleon_google_docstring = False
napoleon_numpy_docstring = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
}

html_theme = "furo"
html_title = f"pysanity {release}"
exclude_patterns = ["_build"]
