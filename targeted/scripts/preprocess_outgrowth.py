# %%

import altair as alt
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted")
hp.setup_altair()

project.log.info("BEGIN: preprocessing barcode outgrowth data.")


# %%
df = hp.data.load_pycashier_outs(
    project.paths.data / "counts-outgrowth.tsv"
).with_columns(pl.col("sample").str.replace("1KB3-P", "").alias("passage"))
library = hp.data.load_pycashier_outs(project.paths.data / "counts-library.tsv")
# %%
df = df.with_columns(pl.col("barcode").is_in(library["barcode"]).alias("in library"))
project.log.info(df.head())


# %%
@project.save_plot()
def outgrowth_barcodes_in_library(df):
    data = df.group_by("passage", "in library").agg(
        pl.col("percent").sum(), pl.col("barcode").n_unique().alias("n_barcodes")
    )

    base = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            alt.X("passage").axis(ticks=False, labelAngle=0),
            color=alt.Color("in library"),
        )
        .properties(height=150)
    )

    percent_bar = base.encode(alt.Y("percent").scale(domain=(0, 100)))

    n_barcode_bar = base.encode(alt.Y("n_barcodes").title("unique barcodes"))
    return percent_bar | n_barcode_bar


outgrowth_barcodes_in_library(df)
# %%
project.log.info(f"writing data to {project.paths.data / 'outgrowth.csv'}")

(
    df.filter(pl.col("in library"))
    .drop("in library")
    .with_columns(
        (pl.col("count") / pl.col("count").sum().over("sample") * 100).alias("percent")
    )
    .write_csv(project.paths.data / "outgrowth.csv")
)
# %%
project.log.info("END: preprocessing barcode treated data.")
# %%
