import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[3]
CLI = REPO / "bin" / "Sanity"
COUNT_TABLE = REPO / "tests" / "count_table.tsv"


@pytest.fixture(scope="session")
def counts() -> pd.DataFrame:
    """Upstream test matrix, genes x cells, including the 17 genes without counts."""
    return pd.read_csv(COUNT_TABLE, sep="\t", index_col=0)


@pytest.fixture(scope="session")
def cli_output(tmp_path_factory):
    """Run the CLI once per method; float64 npy for delta/d_delta/mu/variance, extended text for d_mu."""
    if not CLI.exists():
        pytest.skip(f"{CLI} not found; build it with `mise run build-cli`")
    cache: dict[str, dict[str, Any]] = {}

    def run(method: str) -> dict[str, Any]:
        if method in cache:
            return cache[method]
        npy_dir = tmp_path_factory.mktemp(f"cli_npy_{method}")
        txt_dir = tmp_path_factory.mktemp(f"cli_txt_{method}")
        base = [str(CLI), "-f", str(COUNT_TABLE), "-v_m", method, "-n", "2"]
        subprocess.run(
            [*base, "-d", f"{npy_dir}/", "--npy"], check=True, capture_output=True
        )
        subprocess.run(
            [*base, "-d", f"{txt_dir}/", "-e", "1"], check=True, capture_output=True
        )
        cache[method] = {
            "gene_ids": (npy_dir / "geneID.txt").read_text().split(),
            "log_fold_change": np.load(npy_dir / "delta.npy"),
            "log_fold_change_error": np.load(npy_dir / "d_delta.npy"),
            "mean_ltq": np.load(npy_dir / "mu.npy").ravel(),
            "gene_variance": np.load(npy_dir / "variance.npy").ravel(),
            "mean_ltq_error": np.loadtxt(txt_dir / "d_mu.txt"),
            "variance_likelihood": pd.read_csv(
                txt_dir / "likelihood.txt", sep="\t", index_col=0
            ).to_numpy(),
        }
        return cache[method]

    return run
