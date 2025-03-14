# %%
import altair as alt
import helpers as hp
import polars as pl
import scanpy as sc

# %%
project = hp.Project("scrna/downstream", "scrna.downstream.misc", outs="outs/misc")
project.setup_scanpy()
hp.setup_altair()
dataset = hp.data.get_sc_dataset("231-1KB3")
# %%
adata = sc.read(dataset.final)
# %%
adata = adata[adata.obs["clone-cluster"].isin(["1", "2"])].copy()
adata.obs["sample/clone-cluster"] = (
    adata.obs["sample"].astype("str") + "/" + adata.obs["clone-cluster"].astype("str")
).astype("category")


# %%
sc.pl.correlation_matrix(
    adata, groupby="sample/clone-cluster", save="_sample-clone-cluster.svg"
)


# %%
@project.save_plot()
def post_barcode_totals_by_sample_dotplot(obs=hp.get_obs(adata)):
    post_barcodes = obs.filter(
        pl.col("group").is_in(["250", "400"]) & (pl.col("n_barcodes") == 1)
    )["bgID"].unique()

    totals = (
        obs.filter(pl.col("bgID").is_in(post_barcodes))
        .group_by("bgID")
        .agg(pl.len().alias("n_cells"))
    )
    sample_totals = (
        obs.filter(pl.col("bgID").is_in(post_barcodes))
        .group_by("bgID", "sample")
        .agg(pl.len().alias("n_cells"))
        .with_columns(
            ((pl.col("n_cells") / pl.col("n_cells").sum().over("sample")) * 100).alias(
                "percent"
            )
        )
    )

    # TODO: sample?
    barcode_order = totals.sort("n_cells", descending=True)["bgID"]

    return (
        alt.Chart(sample_totals)
        .mark_point(stroke="black")
        .encode(
            alt.X("sample").axis(labelAngle=300, offset=5),
            alt.Y("bgID").axis(offset=5).sort(barcode_order).title("barcode"),
            alt.Size("percent"),
        )
    )


post_barcode_totals_by_sample_dotplot()
# %%
# ---
