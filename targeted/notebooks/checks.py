# %%
import altair as alt
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted")
hp.setup_altair()


# %%


def load_data(f):
    split_col = pl.col("sample").str.split("-").list
    df = (
        pl.read_csv(f, separator="\t")
        .with_columns(
            split_col.get(0).alias("job"),
            split_col.slice(1).list.join("-").alias("name"),
        )
        .sort("job", "name")
        .select("barcode", "sample", "count", "percent", "job", "name")
    )
    return df


df = pl.concat(
    [
        load_data(project.paths.data / f"counts-{name}.tsv")
        for name in ("library", "treated", "esam", "outgrowth")
    ]
)

df.head()

# %%
library_barcodes = df.filter(pl.col("name") == "libB3")["barcode"].unique()
df = df.with_columns(pl.col("barcode").is_in(library_barcodes).alias("in_library"))


# %%


def plot_count(df, job):
    data = (
        df.filter(pl.col("job") == job)
        .group_by("name", "job", "in_library")
        .agg(pl.col("count").sum())
    )
    x = alt.X("count").scale(domain=[0, 10e6])
    if job != "JA23388":
        x = x.axis(None)
    return (
        alt.Chart(
            data,
            title=alt.Title(job, fontSize=10, orient="right", angle=0),
            height=15 * data["name"].n_unique(),
        )
        .mark_bar()
        .encode(
            alt.Y("name").axis(ticks=False, domain=False).title(None),
            x,
            alt.Color("in_library"),
        )
    )


alt.vconcat(*(plot_count(df, j) for j in sorted(df["job"].unique())), spacing=10)

# %%


def plot_barcodes(df, job):
    data = (
        df.filter(pl.col("job") == job)
        .group_by("name", "job", "in_library")
        .agg(pl.col("barcode").n_unique())
    )
    x = alt.X("barcode").title("unique barcodes")
    if job != "JA23388":
        x = x.axis(None)
    return (
        alt.Chart(
            data,
            title=alt.Title(job, fontSize=10, orient="right", angle=0),
            height=15 * data["name"].n_unique(),
        )
        .mark_bar()
        .encode(
            alt.Y("name").axis(ticks=False, domain=False).title(None),
            x,
            alt.Color("in_library"),
        )
    )


alt.vconcat(*(plot_barcodes(df, j) for j in sorted(df["job"].unique())), spacing=10)
