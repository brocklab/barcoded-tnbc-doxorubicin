# %%
import altair as alt
import helpers as hp
import polars as pl
import scanpy as sc

# %%
project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
project.setup_scanpy()
hp.setup_altair()
dataset = hp.data.get_sc_dataset("231-1KB3")
# %%
adata = sc.read(dataset.final)
# %%
obs = hp.get_obs(adata)
# %%

# %%

# in more than 15 cells in at least two recovered samples
divergent_barcodes = (
    (
        obs.filter(pl.col("group").is_in(("250", "400")) & (pl.col("n_barcodes") == 1))
        .group_by("barcode", "sample")
        .len()
        .rename({"len": "n_cells"})
        .filter(pl.col("n_cells") > 15)
    )
    .group_by("barcode")
    .agg(pl.col("sample").n_unique())
    .filter(pl.col("sample") > 1)["barcode"]
)

# %%


# %%
# %%
@project.save_plot()
def divergent_barcodes_n_cells_bar(obs, divergent_barcodes):
    source = (
        obs.filter(pl.col("barcode").is_in(divergent_barcodes))
        .group_by("barcode", "sample")
        .len()
        .rename({"len": "n_cells"})
    )
    color_scale = hp.plotting.get_alt_palette(
        "sample", groups=source["sample"].unique().to_list()
    )
    return (
        alt.Chart(source)
        .mark_bar()
        .encode(
            alt.X("n_cells").title("cells (N)"),
            alt.Y("barcode"),
            alt.Color("sample").scale(**color_scale),
        )
        .configure_axisY(domain=False, labelFont="monospace")
        .properties(height=100, width=150)
    )


divergent_barcodes_n_cells_bar(obs, divergent_barcodes)
# %%
