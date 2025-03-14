import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

import helpers as hp

project = hp.Project("scrna/downstream", outs="outs/qc")


def add_var_categories(adata):
    adata.var_names_make_unique()
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))


def plot_qc_plots(adata, name, suffix=""):
    hist = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    hist.savefig(project.paths.outs / f"distplot_counts_{name}{suffix}.svg")

    sc.pl.violin(adata, "pct_counts_mt", save=f"_pct_counts_mt_{name}{suffix}.svg")
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        save=f"_counts_{name}{suffix}.svg",
    )


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
