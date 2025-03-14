import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

import helpers as hp

project = hp.Project("scrna/downstream", "scrna.normalize", outs="outs/normalize")


def normalize_counts(adata):
    """normalize counts with shifted logarithm"""
    project.log.info("normalizing counts")

    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)


def plot_normalized_counts(adata, dataset):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm")

    fig.savefig(project.paths.outs / f"normalized_counts_{dataset}.svg")
