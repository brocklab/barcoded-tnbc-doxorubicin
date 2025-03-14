# %%
import sys

import altair as alt
import helpers as hp
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scanpy as sc
import seaborn as sns

alt.data_transformers.enable("vegafusion")

# %%
project = hp.Project(
    "scrna/downstream",
    "scrna.downstream.barcodes",
    outs="outs/deg-survivorship",
)
project.setup_scanpy()
hp.setup_altair()
# %%
adata = sc.read(hp.data.get_sc_dataset("231-1KB3").final)
adata = adata[adata.obs["clone-specific"]].copy()


# %%
def load_degs(sample, cluster, random=False):
    rand = "-randomized" if random else ""
    return (
        pl.read_csv(
            project.paths.data
            / "degs"
            / f"deg-{sample}-cluster-{cluster}-survivorship-low-vs-high{rand}.csv"
        )
    ).with_columns(
        cluster=pl.lit(cluster), sample=pl.lit(sample), randomized=pl.lit(random)
    )


degs = pl.concat(
    load_degs(sample, cluster, random)
    for cluster in ("1", "2")
    for sample in ("Pretreatment", "Dox-550-30HR")
    for random in (True, False)
)


# %%
def volcano_plot(dedf, pval_cutoff=0.05, logFC_cutoff=0.25):
    if dedf["pval"].min() == 0:
        dedf = dedf.with_columns(pl.col("pval") + sys.float_info.min)
    fig, ax = plt.subplots()

    line_kwargs = dict(color="black", linestyle="--")

    sig_genes = pl.col("pval") < pval_cutoff
    up_genes = (pl.col("log2fc") > logFC_cutoff) & sig_genes
    down_genes = (pl.col("log2fc") < -1 * logFC_cutoff) & sig_genes
    else_genes = ~(up_genes | down_genes)

    for f, color in ((up_genes, "blue"), (down_genes, "red"), (else_genes, "grey")):
        ax.scatter(
            x=dedf.filter(f)["log2fc"],
            y=dedf.filter(f)["pval"].log10() * -1,
            c=color,
            s=5,
            alpha=0.8,
        )
    ax.annotate(
        f"{dedf.filter(up_genes).shape[0]} up-regulated",
        (0.95, 1.01),
        xycoords="axes fraction",
        ha="right",
    )
    ax.annotate(
        f"{dedf.filter(down_genes).shape[0]} down-regulated",
        (0.05, 1.01),
        xycoords="axes fraction",
        ha="left",
    )

    x_lim = dedf["log2fc"].abs().max() * 1.1
    ax.set_xlim(-1 * x_lim, x_lim)
    ax.axhline(-1 * np.log10(pval_cutoff), **line_kwargs)
    ax.axvline(logFC_cutoff, **line_kwargs)
    ax.axvline(-1 * logFC_cutoff, **line_kwargs)
    ax.set(xlabel=r"$\log_{2}(FC)$", ylabel=r"$-\log_{10}(pval)$")
    return fig


#


# %%
def save_volcano_plot(degs, cluster, sample, random):
    volcano_plot(
        degs.filter(
            (pl.col("cluster") == cluster),
            (pl.col("sample") == sample),
            (pl.col("randomized") == random),
        )
    )
    rand = "-randomized" if random else ""
    plt.savefig(
        project.paths.outs
        / f"pt-survivorship-{sample}-cluster-{cluster}-deg{rand}-volcano.png",
        dpi=300,
    )


for cluster in "1", "2":
    for sample in "Pretreatment", "Dox-550-30HR":
        for random in True, False:
            save_volcano_plot(degs, cluster, sample, random)


# %%


def get_top_genes(degs, sample, cluster, n_genes=5):
    return (
        degs.filter(
            pl.col("cluster") == cluster,
            pl.col("sample") == sample,
            # ignore ensemble ID only genes
            ~pl.col("gene").str.starts_with("ENSG"),
            ~pl.col("randomized"),
        )
        .sort("pval")
        .head(n_genes)["gene"]
        .to_list()
    )


def get_violin_df(
    adata,
    degs,
    sample,
    cluster,
):
    genes = get_top_genes(degs, sample, cluster)
    print(genes)
    id_cols = ["survivorship", "sample", "clone-cluster"]
    obs_df = (
        sc.get.obs_df(
            adata,
            keys=[
                "survivorship",
                "sample",
                "clone-cluster",
                *genes,
            ],
        )
        .reset_index()
        .rename(columns={"index": "cell"})
    )
    return (
        pl.from_pandas(obs_df)
        .filter(pl.col("sample") == sample, pl.col("clone-cluster") == cluster)
        .unpivot(
            index=["cell", *id_cols], variable_name="gene", value_name="expression"
        )
    )


def plot_top_genes_violinplot(adata, degs, sample, cluster):
    # make sure matplotlib uses a new plot
    plt.figure()

    ax = sns.violinplot(
        data=get_violin_df(adata, degs, sample, cluster),
        x="gene",
        y="expression",
        hue="survivorship",
        palette=hp.plotting.get_palette("survivorship"),
        split=True,
        gap=0.1,
        inner="quart",
        density_norm="count",
    )
    # make it just tall enough to fit the legend on all versions
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    ax.set_ylim(y_min, y_max + 0.05 * y_range)

    ax.set_title(f"{sample}, clone-cluster {cluster} top DEGs")

    ax.legend(title="survivorship", frameon=False, loc="upper left")

    plt.savefig(project.paths.outs / f"deg-top-genes-violin-{sample}-{cluster}.svg")


plot_top_genes_violinplot(adata, degs, "Pretreatment", "2")
# %%
for cluster in "1", "2":
    for sample in "Pretreatment", "Dox-550-30HR":
        plot_top_genes_violinplot(adata, degs, sample, cluster)


# %%
# ---
