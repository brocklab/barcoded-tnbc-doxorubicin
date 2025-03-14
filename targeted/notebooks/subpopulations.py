# TODO: cleanup top - to - bottom
# %%
import altair as alt
import helpers as hp
import polars as pl

# %%

project = hp.Project("targeted", "targeted.treated")
hp.setup_altair()
sample_order = [
    "PARENT-0712",
    "PARENT-0719",
    "PARENT-0729",
    "ESAM-NEG-0729",
    "ESAM-POS-0729",
    "ESAM-NEG-0801",
    "ESAM-POS-0801",
    "ESAM-NEG-2-0809",
    "ESAM-POS-2-0809",
    "PARENT-0810",
]
# %%
df = pl.read_csv(project.paths.data / "counts-esam.tsv", separator="\t")

# fix up sample names
df = df.with_columns(
    pl.col("sample")
    .str.replace("JA23388-", "")
    .str.replace("1KB3-810", "1KB3-0810")
    .str.replace_many({"1KB3-EN": "ESAM-NEG-", "1KB3-EP": "ESAM-POS-"})
    .str.replace("1KB3", "PARENT")
    .str.replace("--", "-")
)

obs = hp.data.obs_final.load()

# lib = pl.read_csv(project.paths.data / "allowlist.csv")
bgIDs = pl.read_csv(project.paths.data / "barcode-bgID.csv")
df = (
    df.join(bgIDs, on="barcode")
    .select("barcode", "bgID", "count", "sample")
    .join(
        obs.filter(pl.col("clone-specific"))
        .select("barcode", "clone-cluster")
        .unique(),
        on="barcode",
    )
    .with_columns(
        (pl.col("count") / pl.col("count").sum().over("sample") * 100).alias("percent")
    )
)


# %%
@project.save_plot()
def esam_separated_totals(df):
    chart = alt.Chart().mark_bar(stroke="black", color="white")
    y = alt.Y("sample").sort(sample_order)
    total = chart.encode(
        alt.X("count").title("total reads"), y.axis(ticks=False, domain=False)
    ).properties(data=df.group_by("sample").agg(pl.col("count").sum()))
    barcodes = chart.encode(
        y.axis(None),
        alt.X("n_barcodes").title("unique barcodes"),
    ).properties(
        data=df.group_by("sample").agg(pl.col("barcode").n_unique().alias("n_barcodes"))
    )
    return (total | barcodes).configure_axis(offset=5)


esam_separated_totals(df)


# %%
@project.save_plot()
def clone_cluster_per_sample_subpopulations(df_esam):
    return (
        alt.Chart(
            df_esam.group_by("sample", "clone-cluster").agg(pl.col("percent").sum())
        )
        .mark_bar(cornerRadius=2)
        .encode(
            alt.X("percent").scale(domain=[0, 100]),
            alt.Y("sample").sort(sample_order).axis(ticks=False, domain=False),
            alt.Color("clone-cluster:N").scale(
                **hp.plotting.get_alt_palette("clone-cluster")
            ),
        )
        .configure_axis(offset=5)
    )


clone_cluster_per_sample_subpopulations(df)

# %%


essential = (
    df.filter(
        pl.col("sample").is_in(
            [
                "PARENT-0712",
                # "PARENT-0719",
                # "PARENT-0729",
                "ESAM-NEG-0729",
                "ESAM-POS-0729",
                # "ESAM-NEG-0801",
                # "ESAM-POS-0801",
                "ESAM-NEG-2-0809",
                "ESAM-POS-2-0809",
                # "PARENT-0810",
            ]
        )
    )
).with_columns(
    pl.col("sample")
    .replace_strict(
        old=[
            "PARENT-0712",
            "ESAM-NEG-0729",
            "ESAM-POS-0729",
            "ESAM-NEG-2-0809",
            "ESAM-POS-2-0809",
        ],
        new=["0", "1", "1", "2", "2"],
    )
    .alias("sorts"),
    pl.col("sample")
    .replace_strict(
        old=[
            "PARENT-0712",
            "ESAM-NEG-0729",
            "ESAM-POS-0729",
            "ESAM-NEG-2-0809",
            "ESAM-POS-2-0809",
        ],
        # new=["library", "ESAM-NEG", "ESAM-POS", "ESAM-NEG", "ESAM-POS"],
        new=["library", "ESAM-", "ESAM+", "ESAM-", "ESAM+"],
    )
    .alias("sample"),
)


# %%
@project.save_plot()
def frequency_clone_cluster_sorted_bar_dot(essential):
    shared_encodings = (
        alt.X("sample").title(None).axis(offset=20),
        alt.Color("clone-cluster").scale(
            **hp.plotting.get_alt_palette("clone-cluster")
        ),
    )

    sort_col = (
        alt.Column("sorts")
        .spacing(5)
        .title("sorts (N)")
        .header(
            orient="bottom",
            titlePadding=-40,
            titleOrient="right",
            titleAnchor="start",
            titleBaseline="line-bottom",
            titleLineHeight=100,
            titleAngle=0,
            labelPadding=-55,
        )
    )
    bar = (
        alt.Chart(essential)
        .mark_bar()
        .encode(
            *shared_encodings,
            alt.Y("percent").scale(domain=[0, 100]).axis(offset=5),
            sort_col.title(None),
        )
    )
    dots = (
        alt.Chart(
            essential.group_by("sample", "clone-cluster", "sorts").agg(
                pl.col("bgID").n_unique().alias("unique barcodes")
            )
        )
        .mark_circle(size=50)
        .encode(*shared_encodings, alt.Y("unique barcodes").axis(offset=5), sort_col)
    )

    return (
        alt.hconcat(
            *[
                c.resolve_scale(x="independent").properties(
                    spacing=4, height=200, width=alt.Step(25)
                )
                for c in (bar, dots)
            ],
            spacing=5,
        )
        .configure_legend(orient="none", legendX=340, legendY=0)
        .configure_view(strokeWidth=1, stroke="black")
    )


frequency_clone_cluster_sorted_bar_dot(essential)
# %%
