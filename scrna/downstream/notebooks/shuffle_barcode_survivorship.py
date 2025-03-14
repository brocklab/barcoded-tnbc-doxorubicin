# %%
import helpers as hp
import polars as pl
import scanpy as sc

# %%
project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
project.setup_scanpy()
dataset = hp.data.get_sc_dataset("231-1KB3")
# %%
adata = sc.read(dataset.final)
obs = hp.get_obs(adata)
# %%
randomized_survivorship = (
    obs.filter((pl.col("n_barcodes") == 1) & (pl.col("survivorship").is_not_null()))
    .select("barcode", "survivorship")
    .unique()
    # sort needed to get consistent results
    .sort("barcode")
    .with_columns(
        pl.col("survivorship").shuffle(seed=231).alias("randomized-survivorship")
    )
    .select("barcode", "randomized-survivorship")
)
# %%
# %%
randomized_survivorship.write_csv(
    project.paths.data / "randomized-barcode-survivorship.csv"
)
