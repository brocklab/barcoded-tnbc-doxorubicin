# TODO: delete
# %%
import logging

import altair as alt
import helpers as hp
import matplotlib.pyplot as plt
import polars as pl
import scanpy as sc

# %%
hp.set_project("scrna/downstream")
hp.setup_altair()
hp.setup_scanpy()
logger = logging.getLogger("scrna.downstream.figures")
# %%
samples = {
    "PT": "Pretreatment",
    "Dox": "Dox-550-30h",
    "A4": "Dox-250-4",
    "A6": "Dox-250-6",
}
dataset = hp.data.get_sc_dataset("231-1KB3")
adata = sc.read(dataset.barcoded)
obs = hp.get_obs(adata).with_columns(
    pl.col("sample").cast(str).replace(samples), pl.col("sample").alias("sample-id")
)


# %%
@hp.save_plot()
def cell_cycle_by_sample(obs):
    data = (
        obs.group_by("sample", "phase")
        .agg(pl.len().alias("n_cells"))
        .with_columns(
            (pl.col("n_cells") / pl.col("n_cells").sum())
            .over("sample")
            .alias("proportion")
        )
    )

    return (
        alt.Chart(data)
        .mark_bar()
        .encode(
            alt.X("proportion").axis(offset=5).title("proportion"),
            alt.Y("sample")
            .sort(list(samples.values()))
            .axis(domain=False, ticks=False, offset=5),
            color=alt.Color("phase:N"),
        )
        .properties(height=150)
    )


cell_cycle_by_sample(obs)


# %%
def plot_n_barcodes(adata):
    adata = adata.copy()
    adata.obs["n_barcodes"] = adata.obs["n_barcodes"].fillna(0)
    og = adata.obs["n_barcodes"]
    fig, axes = plt.subplots(2, 2)

    for i, ax in enumerate(axes.flat):
        adata.obs["n_barcodes"] = og
        adata.obs.loc[adata.obs["n_barcodes"] != i, "n_barcodes"] = float("nan")
        adata.obs["n_barcodes"] = adata.obs["n_barcodes"].astype("category")
        sc.pl.umap(
            adata,
            color="n_barcodes",
            frameon=False,
            legend_loc=None,
            title=f"barcodes/cell = {i}",
            ax=ax,
            show=False,
        )
    fig.savefig(hp.paths.outs / "n_barcodes_per_cell.png")
    return fig


plot_n_barcodes(adata)


# %%
def plot_cell_cycle(df):
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            alt.Y("sample").axis(ticks=False),
            alt.X("len").title("# of cells"),
            facet=alt.Facet("barcode").header(labelFont="monospace").columns(4),
            color=alt.Color("phase:N"),
        )
        .resolve_scale(
            x="independent",
        )
    ).properties(width=100)


post_barcodes = (
    obs.filter(
        (pl.col("sample").is_in(("Dox-250-4", "Dox-250-6")))
        & (pl.col("barcode") != "unknown")
    )["barcode"]
    .unique()
    .to_list()
)

df_cell_cycle = (
    obs.filter(pl.col("barcode").is_in(post_barcodes))
    .group_by("sample", "barcode", "phase", "clone-cluster")
    .agg(pl.len())
)
# %%
(p := plot_cell_cycle(df_cell_cycle)).save(
    hp.paths.outs / "cell-cycle-totals-250-barcodes.svg"
)
p
# %%
(
    p := plot_cell_cycle(
        df_cell_cycle.filter(pl.col("sample").is_in(["Pretreatment", "Dox-550-30h"]))
    )
).save(hp.paths.outs / "cell-cycle-totals-250-barcodes-without-250.svg")
p
# %%
