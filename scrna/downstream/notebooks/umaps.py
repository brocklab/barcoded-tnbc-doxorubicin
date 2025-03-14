# %%
import helpers as hp
import matplotlib.pyplot as plt
import polars as pl
import scanpy as sc
from helpers.plotting.umap import Umap

# %%
project = hp.Project("scrna/downstream", outs="outs/umaps")
project.setup_scanpy()
# %%
adata = sc.read(hp.data.get_sc_dataset("231-1KB3").final)

# %%
Umap(adata, "sample").plot(
    legend_kwargs=dict(
        bbox_to_anchor=(0.95, 0.5),
        loc="center left",
    ),
    use_uns_colors=True,
    save=project.paths.outs / "umap-sample.png",
)
# %%
# umaps of sample with each group added one at time


def umap_by_groups(adata_, groups):
    adata_ = adata_.copy()
    adata_.obs.loc[~adata_.obs["group"].isin(groups), "sample"] = float("nan")
    group_str = "-".join(groups)
    Umap(adata_, "sample").plot(
        legend_kwargs=dict(
            bbox_to_anchor=(0.95, 0.5),
            loc="center left",
        ),
        use_uns_colors=True,
        save=project.paths.outs / f"umap-sample-grps-{group_str}.png",
    )


umap_by_groups(adata, ["PT"])
# %%
grps = []
for group in ["PT", "30HR", "250", "400"]:
    grps.append(group)
    umap_by_groups(adata, grps)


# %%
Umap(adata, "group").plot(
    legend_kwargs=dict(
        bbox_to_anchor=(0.95, 0.5),
        loc="center left",
    ),
    use_uns_colors=True,
    save=project.paths.outs / "umap-group.png",
)
# %%
Umap(adata, "leiden").plot(
    legend_kwargs=dict(
        bbox_to_anchor=(0.95, 0.5),
        loc="center left",
    ),
    save=project.paths.outs / "umap-leiden.png",
)
# %%
Umap(adata, "ESAM").plot(save=project.paths.outs / "umap-esam.png")
# %%
(
    Umap(adata, "survivor").plot(
        cmap=dict(yes="blue", no="orange"),
        save=project.paths.outs / "umap-survivor.png",
    )
)
# %%
Umap(adata, "clone-cluster").plot(
    cmap=hp.plotting.get_palette("clone-cluster"),
    title="clone cluster",
    save=project.paths.outs / "umap-clone-cluster.png",
)
# %%
Umap(adata, "phase").plot(save=project.paths.outs / "umap-phase.png")
# %%
Umap(adata, "total_counts").plot(save=project.paths.outs / "umap-total-counts.png")
# %%
# individual umaps for each clone cluster in pretreatment and Dox-550-30HR
sample_filter = adata.obs["sample"].isin(["Pretreatment", "Dox-550-30HR"])
two_sample_legend = dict(
    bbox_to_anchor=(0.04, 0.90),
    loc="center left",
)
for cluster in "1", "2":
    adata_umap = adata.copy()
    adata_umap.obs.loc[
        ~sample_filter | ~(adata.obs["clone-cluster"] == cluster), "sample"
    ] = float("nan")
    (
        Umap(adata_umap, key="sample").plot(
            title=f"clone-cluster {cluster}",
            save=project.paths.outs / f"umap-clone-cluster-{cluster}-pt-dox.png",
            legend_kwargs=two_sample_legend,
        )
    )
# %%
# individual umaps for each clone cluster in pretreatment and Dox-550-30HR
sample_filter = adata.obs["sample"].isin(["Pretreatment"])
two_sample_legend = dict(
    bbox_to_anchor=(0.04, 0.90),
    loc="center left",
)
obs = "survivorship"
for cluster in "1", "2":
    adata_umap = adata.copy()
    cell_filters = ~sample_filter | ~(adata.obs["clone-cluster"] == cluster)
    adata_umap.obs.loc[cell_filters, obs] = float("nan")
    adata_umap.obs.loc[adata_umap.obs["barcode"] == "unknown", obs] = float("nan")
    adata_umap.obs[obs] = adata_umap.obs[obs].cat.set_categories(["high", "low"])
    (
        Umap(adata_umap, key=obs).plot(
            cmap=hp.plotting.get_palette("survivorship", ["high", "low"]),
            title=f"survivorship clone-cluster {cluster}",
            save=project.paths.outs / f"umap-survivorship-clone-cluster-{cluster}.png",
            legend_kwargs=two_sample_legend,
        )
    )

# %%
adata_umap = adata.copy()
sample_filter = adata.obs["sample"].isin(["Pretreatment"])
two_sample_legend = dict(
    bbox_to_anchor=(0.04, 0.90),
    loc="center left",
)
obs = "survivorship"
adata_umap = adata.copy()
cell_filters = ~sample_filter
adata_umap.obs.loc[cell_filters, obs] = float("nan")
adata_umap.obs.loc[adata_umap.obs["barcode"] == "unknown", obs] = float("nan")
adata_umap.obs[obs] = adata_umap.obs[obs].cat.set_categories(["high", "low"])
(
    Umap(adata_umap, key=obs).plot(
        cmap=hp.plotting.get_palette("survivorship", ["high", "low"]),
        title="Pretreatment survivorship",
        save=project.paths.outs / "umap-survivorship-Pretreatment.png",
        legend_kwargs=two_sample_legend,
    )
)


# %%
samples = adata.obs["sample"].cat.categories.to_list()
adata_plot = adata.copy()
fig, axes = plt.subplots(2, 4, figsize=(8, 4))
for sample, ax in zip(samples, axes.flat):
    adata_plot.obs["sample"] = adata.obs["sample"]
    adata_plot.obs.loc[adata.obs["sample"] != sample, "sample"] = float("nan")
    sc.pl.umap(
        adata_plot,
        color="sample",
        title=sample,
        legend_loc=None,
        frameon=False,
        show=False,
        ax=ax,
    )

fig.savefig(project.paths.outs / "umap-per-sample.png", bbox_inches="tight", dpi=300)
fig
# %%
barcodes = (
    hp.get_obs(adata)
    .filter(~pl.col("barcode").is_null())
    .group_by("barcode")
    .agg(pl.len())
    .sort("len", descending=True)["barcode"][:10]
    .to_list()
)
adata_plot = adata.copy()
adata_plot.obs["barcode"] = adata.obs["barcode"].cat.set_categories(barcodes)
Umap(adata_plot, key="barcode").plot(
    legend=False,
    title="top 10 barcodes",
    save=project.paths.outs / "umap-top-10-barcodes.png",
)
# %%
hp.get_obs(adata).group_by("barcode", "clone-cluster").agg(pl.len()).sort(
    "len", descending=True
)
# %%
adata_plot = adata.copy()
adata_plot.obs.loc[~adata.obs["group"].isin(["PT", "30HR"]), "sample"] = float("nan")
Umap(adata_plot, "sample").plot(
    # legend_kwargs=dict(
    #     bbox_to_anchor=(0.95, 0.5),
    #     loc="center left",
    # ),
    use_uns_colors=True,
    # title=barcode,
    save=project.paths.outs / "umap-sample-pt-dox-30hr-only.png",
)
# %%
