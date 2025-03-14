# %%
import altair as alt
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted")
hp.setup_altair()

# %%
df = pl.read_csv(project.paths.data / "outgrowth.csv")
project.log.info(df.head())


# %%
@project.save_plot()
def unique_barcodes_outgrowth(df):
    return (
        alt.Chart(
            df.group_by("passage").agg(pl.col("barcode").n_unique().alias("n_barcodes"))
        )
        .mark_bar()
        .encode(
            alt.X("passage:O").axis(ticks=False, labelAngle=0),
            alt.Y("n_barcodes").title("unique barcodes"),
        )
    )


unique_barcodes_outgrowth(df)


# %%
@project.save_plot()
def top_ten_P31_barcodes_abundance_change(df):
    top_ten_P31_barcodes = df.filter(pl.col("passage") == 31).sort(
        "percent", descending=True
    )["barcode"][:10]

    base = (
        alt.Chart(df.filter(pl.col("barcode").is_in(top_ten_P31_barcodes)))
        .encode(
            alt.X("passage:O"),
            alt.Y("percent"),
            color=alt.Color("barcode").legend(labelFont="monospace"),
        )
        .properties(height=200, width=200)
    )

    dots = base.mark_circle()
    line = base.mark_line()

    return dots + line


top_ten_P31_barcodes_abundance_change(df)
# %%
