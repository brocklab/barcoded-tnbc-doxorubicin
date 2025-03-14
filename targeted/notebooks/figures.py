# %%
import altair as alt
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted", "targeted.treated")
hp.setup_altair()
# %%
project.log.info("BEGIN: analyzing treated targeted barcodes.")

# %%
df = hp.data.targeted_final.load()

sample_order = [
    "PT-0",
    *sorted(
        (
            *df.filter(pl.col("sample") != "PT-0")["sample"].unique().sort().to_list(),
            "550-4",
        )
    ),
]
# the ranking doesn't special case 'unknown' in any way
df = df.with_columns(
    pl.col("percent").rank(descending=True).over(pl.col("sample")).alias("rank")
)

post_treatment_barcodes = df.filter(
    (pl.col("group") != "PT") & (pl.col("barcode") != "unknown")
)["barcode"].unique()

obs = hp.data.obs_final.load()

clone_cluster_to_esam = {"1": "ESAM+", "2": "ESAM-"}
barcode_esam = (
    obs.filter(pl.col("sample") == "Pretreatment")
    .select("barcode", "clone-cluster")
    .unique()
    .with_columns(
        pl.col("clone-cluster").cast(str).replace(clone_cluster_to_esam).alias("esam")
    )
    # for now remove multiple barcodes
    .filter(pl.col("barcode").cast(str).str.len_chars() < 25)
)

growth = (
    pl.read_csv(project.paths.data / "outgrowth.csv")
    .join(barcode_esam, on="barcode", how="left")
    .fill_null("unknown")
)


# %%
@project.save_plot()
def clone_cluster_abundance_passage_area(growth):
    return (
        alt.Chart(growth)
        .mark_area()
        .encode(
            alt.X("passage").axis(values=sorted(growth["passage"].unique().to_list())),
            alt.Detail("barcode"),
            alt.Color("clone-cluster")
            # .scale(**clone_cluster_scale)
            .scale(**hp.plotting.get_alt_palette("clone-cluster"))
            .legend(None),  # will use legend from neighbor,
            alt.Y("percent").scale(domain=[0, 100]),
        )
        .properties(width=200, height=200)
        .configure_axis(offset=10)
    )


@project.save_plot()
def clone_cluster_unique_barcodes_dots(growth):
    source = growth.group_by("clone-cluster", "passage").agg(
        pl.col("barcode").n_unique().alias("unique barcodes")
    )
    return (
        alt.Chart(source)
        .mark_circle(size=100)
        .encode(
            alt.X("passage:O")
            .axis(values=sorted(growth["passage"].unique().to_list()), labelAngle=0)
            .scale(padding=0),
            alt.Color("clone-cluster").scale(
                **hp.plotting.get_alt_palette("clone-cluster")
            ),
            alt.Y("unique barcodes"),
        )
        .properties(width=125, height=200)
        .configure_axis(offset=10)
        .configure_legend(orient="none", legendX=75)
    )


# %%
clone_cluster_abundance_passage_area(growth)
# %%
clone_cluster_unique_barcodes_dots(growth)


# %%
def add_offset(post_n_barcodes, scale=2):
    return post_n_barcodes.with_columns(
        pl.col("n_barcodes").rank("ordinal").over("group", "n_barcodes").alias("rank"),
        pl.col("n_barcodes").count().over("group", "n_barcodes").alias("count"),
    ).with_columns(
        pl.when(pl.col("count") > 1)
        .then(((2 * (pl.col("rank") - 1) / (pl.col("count") - 1)) - 1) / scale)
        .otherwise(0)
        .alias("offset")
    )


@project.save_plot()
def barcodes_post_treatment_per_sample_per_group_swarmplot(df):
    post = (post := df.filter(pl.col("sample") != "PT-0")).join(
        post.group_by("barcode").agg(pl.col("sample").n_unique().alias("n_samples")),
        on="barcode",
    )

    post_n_barcodes = post.group_by(["sample", "group"]).agg(
        pl.col("barcode").n_unique().alias("n_barcodes")
    )
    post_n_barcodes = add_offset(post_n_barcodes, scale=3)
    grp_palette = hp.plotting.get_palette(
        "group", post_n_barcodes["group"].unique().to_list()
    )
    return (
        (
            alt.Chart(post_n_barcodes)
            .mark_circle(size=50, opacity=1)
            .encode(
                alt.X("group")
                .axis(ticks=False, labelAngle=0)
                .title("Doxorubicin (nM)"),
                alt.XOffset("offset").scale(domain=[-1, 1]),
                alt.Y("n_barcodes").title("unique barcodes"),
                color=alt.Color("group")
                .scale(domain=grp_palette.keys(), range=grp_palette.values())
                .legend(None),
            )
        )
        .configure_legend(orient="none", legendX=65, legendY=5)
        .properties(
            # width=15 * post_n_barcodes["sample"].n_unique(),
            width=100,
            height=200,
        )
    )


barcodes_post_treatment_per_sample_per_group_swarmplot(df)

# %%


@project.save_plot()
def percent_unknown(df):
    unknown = pl.concat(
        [
            df.with_columns((pl.col("barcode") == "unknown").alias("unknown")).select(
                "sample",
                "group",
                "percent",
                "unknown",
            ),
            # add back in removed sample C4 (aka C5 now)
            pl.DataFrame(
                dict(sample=["550-5"], group=["550"], percent=[100.0], unknown=[True])
            ),
        ]
    )

    return (
        alt.Chart(unknown.group_by("sample", "unknown").agg(pl.col("percent").sum()))
        .mark_bar(stroke="black", cornerRadius=2)
        .encode(
            alt.X("sample").axis(ticks=False, domain=False).sort(sample_order),
            alt.Y("percent").axis(offset=5).scale(domain=(0, 100)),
            color=alt.Color("unknown").scale(range=["#FFFFFF", "red"]),
        )
        .properties(width=200, height=150)
    )


percent_unknown(df)


# %%


@project.save_plot()
def post_treatment_barcode_dotplot(df):
    group_palette = hp.plotting.get_alt_palette("group", df["group"].unique().to_list())
    group_color_enc = alt.Color("group").scale(**group_palette)

    barcode_order = (
        # pretreatment abundance
        # df.filter(pl.col("sample") == "PT-0").sort("rank")["barcode"].to_list()
        #
        # num of post treatment samples
        df.filter(pl.col("sample") != "PT-0")
        .group_by("barcode")
        .agg(
            pl.col("sample").n_unique(),
            pl.col("group").n_unique(),
            pl.col("percent").sum(),
        )
        # include all three so result is deterministic
        .sort("sample", "group", "percent", "barcode", descending=True)["barcode"]
        .to_list()
    )

    data = df

    dot_dims = dict(width=10 * len(barcode_order), height=150)
    barcode_x_enc = alt.X("barcode:O").sort(barcode_order)
    spacing = 1

    total_samples_bar = (
        alt.Chart(
            data.filter(pl.col("sample") != "PT-0")
            .group_by("barcode", "group")
            .agg(pl.col("sample").n_unique())
        )
        .mark_bar()
        .encode(
            barcode_x_enc.axis(None),
            alt.Y("sample")
            .title(("post-treatment", "samples (N)"))
            .axis(offset=10, titleFontSize=8, tickCount=5),
            group_color_enc,
        )
        .properties(width=dot_dims["width"], height=50)
    )
    dots = (
        alt.Chart(data)
        .mark_circle(xOffset=0.3)
        .encode(
            (
                (barcode_x_enc)
                .axis(
                    domain=False,
                    labels=False,
                    ticks=False,
                    grid=True,
                    gridDash=[3, 1],
                    gridWidth=0.3,
                    gridColor="black",
                    gridOpacity=0.1,
                )
                .scale(type="band", nice=False)
            ),
            (alt.Y("sample").sort(sample_order).axis(offset=10).title(None)),
            alt.Size("percent").scale(rangeMax=150),
            group_color_enc,
        )
        .properties(**dot_dims)
    )
    return alt.vconcat(
        total_samples_bar,
        dots,
        spacing=spacing,
    )


post_treatment_barcodes = df.filter(
    (pl.col("group") != "PT") & (pl.col("barcode") != "unknown")
)["barcode"].unique()

post_treatment_barcode_dotplot(
    df.filter(pl.col("barcode").is_in(post_treatment_barcodes))
)


# %%
@project.save_plot()
def barcode_dotplot(df):
    group_palette = hp.plotting.get_alt_palette("group", df["group"].unique().to_list())
    group_color_enc = alt.Color("group").scale(**group_palette)

    barcode_order = (
        # pretreatment abundance
        df.filter(pl.col("barcode") != "unknown", pl.col("sample") == "PT-0")
        .sort("percent", descending=True)["barcode"]
        .to_list()
    )

    data = df.filter(pl.col("barcode") != "unknown")

    dot_dims = dict(width=1 * len(barcode_order), height=150)
    barcode_x_enc = alt.X("barcode:O").sort(barcode_order)

    dots = (
        alt.Chart(data)
        .mark_circle(xOffset=0.3)
        .encode(
            (
                (barcode_x_enc)
                .axis(
                    domain=False,
                    labels=False,
                    ticks=False,
                )
                .scale(type="band", nice=False)
            ),
            (alt.Y("sample").sort(sample_order).axis(offset=10).title(None)),
            alt.Size("percent").scale(rangeMax=125),
            group_color_enc,
        )
        .properties(**dot_dims)
    )
    return dots


barcode_dotplot(df)


# %%
@project.save_plot()
def top_14_barcode_dot_plot_facet(df):
    top_14 = (
        df.filter((pl.col("group") != "PT") & (pl.col("barcode") != "unknown"))
        .group_by("bgID")
        # sorted by # of post-treatment samples
        .agg(pl.col("sample").n_unique())
        .sort("sample", "bgID", descending=True)[:14]["bgID"]
    )
    source = df.filter(pl.col("bgID").is_in(top_14))
    return (
        alt.Chart(source)
        .mark_circle(size=80)
        .encode(
            alt.X("group").scale(domain=["PT", "250", "400", "550"]).title(None),
            alt.Y("percent").scale(type="log", domain=[0.001, 300]),
            alt.Color("group").scale(**hp.plotting.get_alt_palette("group")),
            alt.Facet("bgID")
            .columns(7)
            .spacing(1)
            .title(None)
            .header(labelPadding=0, labelFontWeight="bold"),
        )
        .properties(height=155, width=80)
        .configure_view(strokeWidth=1, stroke="black")
    )


top_14_barcode_dot_plot_facet(df)


# %% [markdown]
# Is persisting after growth based on initial abundance?
# Or in other words is the abundance slowly skewing from the initial sort?
# %%
@project.save_plot()
def outgrowth_persistence_vs_initial_abundance(growth):
    source = growth.filter(pl.col("passage") == growth["passage"].min()).join(
        growth.group_by("barcode").agg(
            (pl.col("passage").max() - 17).alias("n_passages")
        ),
        on="barcode",
    )
    return (
        alt.Chart(source)
        .mark_circle(opacity=0.5)
        .encode(
            alt.X("percent")
            .scale(type="log")
            .axis(
                labelExpr="datum.value == pow(10, round(log(datum.value) / 2.303)) ? datum.label : ''",
                offset=10,
            ),
            alt.Y("n_passages:O").title("passages (N)").sort("-y"),
            alt.YOffset("jitter:Q"),
        )
        .transform_calculate(
            # Generate Gaussian jitter with a Box-Muller transform
            jitter="sqrt(-2*log(random()))*cos(2*PI*random())"
        )
        .properties(height=200)
    )


outgrowth_persistence_vs_initial_abundance(growth)
# %%
growth.filter(pl.col("passage") == growth["passage"].min()).join(
    growth.group_by("barcode").agg(pl.col("passage").n_unique().alias("n_passages")),
    on="barcode",
)


# %%
@project.save_plot()
def shared_bar_post(df):
    return (
        alt.Chart(df.filter(pl.col("barcode") != "unknown"))
        .mark_bar()
        .encode(
            alt.X("sample"),
            alt.Y("percent")
            .title(["percent of", "total population"])
            .scale(domain=(0, 100)),
            color=alt.Color("n_samples").title("# of samples"),
        )
        .properties(width=300, height=200)
    )


# %%
post = (post := df.filter(pl.col("sample") != "PT-0")).join(
    post.group_by("barcode").agg(pl.col("sample").n_unique().alias("n_samples")),
    on="barcode",
)


shared_bar_post(post)


# %%
@project.save_plot()
def clone_cluster_per_treated_samples(df):
    # add back 550-5
    extra = pl.DataFrame(
        {
            "sample": ["550-5"],
            "percent": [100.0],
            "clone-cluster": ["unknown"],
        }
    )
    source = pl.concat(
        (
            df.join(barcode_esam, on="barcode", how="left")
            .fill_null("unknown")
            .select("sample", "percent", "clone-cluster"),
            extra,
        ),
    )

    return (
        alt.Chart(source)
        .mark_bar()
        .encode(
            alt.Y("sample")
            .axis(domain=False, ticks=False, offset=5)
            .sort(sample_order),
            alt.X("percent").scale(domain=[0, 100]),
            alt.Color("clone-cluster")
            .scale(**hp.plotting.get_alt_palette("clone-cluster"))
            .legend(orient="top", offset=10),
        )
        .properties(width=250, height=225)
    )


clone_cluster_per_treated_samples(df)
# %%
df
# %%
