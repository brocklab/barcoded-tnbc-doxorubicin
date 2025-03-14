# TODO: cleanup
# %%
import helpers as hp
import polars as pl
import scanpy as sc
from helpers.plotting.umap import Umap

# %%
project = hp.Project("scrna/downstream", outs="outs/umaps")
project.setup_scanpy()
hp.setup_altair()
# %%
adata = sc.read(hp.data.get_sc_dataset("231-1KB3").final)

# %%
Umap(adata, "sample").plot(
    legend_kwargs=dict(
        bbox_to_anchor=(0.95, 0.5),
        loc="center left",
    ),
    use_uns_colors=True,
)

# %%
barcode = "TTTCTTCACCGCTATTGAGC"
adata_plot = adata.copy()
# adata_plot.obs.loc[adata.obs["barcode"] != barcode, "sample"] = float("nan")
# %%

# %%
# for sample in hp.get_obs(adata).filter(pl.col("barcode") == barcode)["sample"].unique():
#     print(sample)
#     sc.pl.umap(
#         adata_plot,
#         color="sample",
#         mask_obs=((adata.obs["barcode"] == barcode) & (adata.obs["sample"] == sample)),
#     )

# %%
hp.get_obs(adata).filter(pl.col("barcode") == barcode).group_by("sample").agg(pl.len())
# %%

# adata_plot = adata.copy()
# obs = adata.obs.copy()

# # %%
# barcode = "TTTCTTCACCGCTATTGAGC"
# sample = "Dox-250-4"


# def plot_umap_barcode_sample(adata_, barcode, sample):
#     adata_.obs.loc[
#         ~((adata.obs["barcode"] == barcode) & (adata.obs["sample"] == sample)), "sample"
#     ] = float("nan")
#     sc.pl.umap(
#         adata_plot,
#         color="sample",
#         title=f"{barcode}, {sample}",
#         legend_loc=None,
#         frameon=False,
#     )


# # %%

# for sample in hp.get_obs(adata).filter(pl.col("barcode") == barcode)["sample"].unique():
#     adata_plot.obs = obs.copy()
#     plot_umap_barcode_sample(adata_plot, barcode, sample)
# %%


# %%
barcode_totals = (
    hp.get_obs(adata)
    .filter(pl.col("group").is_in(("250", "400")) & pl.col("barcode").is_not_null())
    .group_by("barcode")
    .agg(pl.len().alias("n_cells"))
    .filter(pl.col("n_cells") > 50)
    .sort("n_cells", descending=True)
)  # .plot.bar(x='n_cells', y='barcode')

# %%
# for now ignore doublet barcodes
post_barcodes = (
    hp.get_obs(adata)
    .filter(pl.col("group").is_in(["250", "400"]) & (pl.col("n_barcodes") == 1))[
        "barcode"
    ]
    .unique()
)

# %%


# %%


# %%

# %%


def plot_umap_barcode_sample(adata_, bgID):
    samples = hp.get_obs(adata_).filter(pl.col("bgID") == bgID)["sample"].unique()
    labels = [f"{bgID}, {sample}" for sample in samples]

    for label, sample in zip(labels, samples):
        adata_.obs[label] = adata_.obs["sample"]
        adata_.obs.loc[
            ~((adata.obs["bgID"] == bgID) & (adata.obs["sample"] == sample)),
            label,
        ] = float("nan")

    # n_cells = ((~adata_.obs['sample'].isna()).sum())
    titles = [f"{label}, n={(~adata_.obs[label].isna()).sum()}" for label in labels]
    sc.pl.umap(
        adata_plot,
        color=labels,
        # title=f"{barcode}, {sample}, n={n_cells}",
        title=titles,
        legend_loc=None,
        frameon=False,
        wspace=0,
        size=20,
        save=f"_{bgID}-samples.png",
    )


adata_plot = adata.copy()
obs = adata.obs.copy()
post_barcodes_ids = (
    hp.get_obs(adata)
    .filter(pl.col("group").is_in(["250", "400"]), pl.col("clone-specific"))
    .group_by("bgID")
    .agg(pl.len())
    .sort("len", descending=True)
    .filter(pl.col("len") >= 50)["bgID"]
)


for bgID in post_barcodes_ids:
    # for sample in hp.get_obs(adata).filter(pl.col('barcode') == barcode)['sample'].unique():
    adata_plot.obs = obs.copy()
    plot_umap_barcode_sample(adata_plot, bgID)


# %%
# superceded by bgIDs
# multi_sample_barcodes = [
#     "GCGAAGGAGCGTTGCGTTGG",
#     "TGTATACCGTCATGTGTGGT",
#     "TTTCTTCACCGCTATTGAGC",
# ]
multi_sample_bgIDs = [
    "bg010",
    "bg032",
    "bg016",
]
# %%


def plot_umap_by_sample_for_barcode(adata_, bgID):
    aplot = adata_.copy()
    aplot.obs.loc[aplot.obs["bgID"] != bgID, "sample"] = float("nan")
    Umap(aplot, key="sample").plot(
        title=bgID,
        save=project.paths.outs / f"umap-{bgID}-by-sample.png",
        legend_kwargs=dict(
            bbox_to_anchor=(0.98, 0.90),
            loc="center right",
        ),
    )


for bgID in multi_sample_bgIDs:
    plot_umap_by_sample_for_barcode(adata, bgID)

# %%


def plot_umap_top_ten_survivors(adata_):
    top_ten_survivor_barcodes = (
        (
            hp.get_obs(adata_)
            .filter(pl.col("n_barcodes") == 1, pl.col("group").is_in(("250", "400")))
            .group_by("barcode")
        )
        .len()
        .sort("len", descending=True)[:10]["barcode"]
        .to_list()
    )
    aplot = adata_.copy()
    aplot.obs.loc[~aplot.obs["barcode"].isin(top_ten_survivor_barcodes), "barcode"] = (
        float("nan")
    )
    Umap(aplot, key="barcode").plot(
        title="top ten survivor clones",
        legend=False,
        save=project.paths.outs / "umap-top-ten-survivor-clones.png",
    )


plot_umap_top_ten_survivors(adata)
# %%
# ---
