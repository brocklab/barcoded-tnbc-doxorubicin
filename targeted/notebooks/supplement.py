# %%
import altair as alt
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted")
hp.setup_altair()

# %%
bfp = pl.read_csv(project.paths.data / "DM2103-bfp.tsv").with_columns(
    pl.when(pl.col("bfp") == 99).then(pl.lit(100)).otherwise(pl.col("bfp")).alias("bfp")
)


# %%
@project.save_plot()
def bfp_positive_sample(bfp):
    return (
        alt.Chart(bfp)
        .mark_bar(stroke="black", color="white")
        .encode(
            alt.Y("sample").axis(ticks=False, offset=5, domain=False),
            alt.X("bfp").title("BFP+ (%)"),
        )
        .properties(height=300)
    )


bfp_positive_sample(bfp)

# %%
