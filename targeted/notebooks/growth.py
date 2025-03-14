# %%
import math

import altair as alt
import helpers as hp
import polars as pl
from scipy.stats import ttest_ind

# %%
project = hp.Project("targeted")
hp.setup_altair()
# %%
treated = hp.data.targeted_final.load()


# %%
def barcode_clone_clusters():
    return (
        hp.data.obs_final.load()
        .filter(pl.col("sample") == "Pretreatment", pl.col("clone-specific"))
        .select("barcode", "clone-cluster")
        .unique()
    )


def load_data(treated):
    initial = treated.filter(pl.col("sample") == "PT-0").select(
        "barcode", "percent", pl.lit(15, pl.Int64).alias("passage")
    )
    rest = pl.read_csv(project.paths.data / "outgrowth.csv").select(
        "barcode", "percent", "passage"
    )

    df = pl.concat((initial, rest))
    df = df.join(barcode_clone_clusters(), on="barcode", how="left").fill_null(
        "unknown"
    )
    df = df.filter(
        pl.col("barcode").is_in(df.filter(pl.col("passage") == 15)["barcode"])
    )
    df = df.with_columns(
        pl.col("barcode")
        .is_in(treated.filter(pl.col("sample") != "PT-0")["barcode"])
        .alias("survivor")
    )
    return df


df = load_data(treated)


# %%
@project.save_plot()
def barcode_dotplot_outgrowth(df):
    color_enc = alt.Color("clone-cluster").scale(
        **hp.plotting.get_alt_palette("clone-cluster")
    )
    barcode_order = (
        # pretreatment abundance
        df.filter(pl.col("passage") == 15)
        .sort("percent", descending=True)["barcode"]
        .to_list()
    )

    data = df.filter(pl.col("barcode") != "unknown")

    dot_dims = dict(height=2.5 * len(barcode_order), width=150)
    barcode_x_enc = alt.Y("barcode:O").sort(barcode_order)

    dots = (
        alt.Chart(data)
        .mark_circle(xOffset=0.3)
        .encode(
            (
                (barcode_x_enc)
                .axis(domain=False, labels=False, ticks=False)
                .scale(type="band", nice=False)
            ),
            alt.X("passage:O").axis(offset=10, labelAngle=0),
            alt.Size("percent").scale(rangeMax=200),
            color_enc,
        )
        .properties(**dot_dims)
    )
    return dots


barcode_dotplot_outgrowth(df)


# %%
def get_growth_rate(initial_percent, final_percent):
    time = 5 * 24 * 2
    N_0 = 500_000
    # assume an average doubling time of 27 hours
    avg_r = math.log(2) / 27
    N_f = N_0 * (math.e ** (time * avg_r))

    # 2 passages, 5 days long each (ignoring possible effects from downsampling)
    # assume abundance after 2 passages is ~ abundance after one passage
    return math.log((N_f * final_percent / 100) / (N_0 * initial_percent / 100)) / time


def get_growth_rates(df, barcode):
    rows = df.filter(pl.col("barcode") == barcode).sort("passage").rows(named=True)
    for i, row in enumerate(rows[:-1]):
        yield get_growth_rate(row["percent"], rows[i + 1]["percent"])


def get_all_growth_rates(df):
    rates = []
    for barcode in (
        df.group_by("barcode")
        .agg(pl.col("passage").n_unique())
        .filter(pl.col("passage") > 2)["barcode"]
    ):
        rates.append(
            pl.DataFrame(get_growth_rates(df, barcode), schema="r").with_columns(
                pl.lit(barcode).alias("barcode")
            )
        )

    return pl.concat(rates)


# %%
rates = (
    get_all_growth_rates(df)
    .group_by("barcode")
    .agg(pl.col("r").mean(), pl.col("r").count().alias("n_rates"))
)

# %%


def post_treatment_samples_per_barcode(treated):
    post = (
        treated.filter(
            pl.col("sample") != "PT-0",
            pl.col("barcode") != "unknown",
        )
        .group_by("barcode")
        .agg(pl.col("sample").n_unique().alias("n-post-samples"))
    )
    pre = (
        treated.filter(
            pl.col("barcode") != "unknown", ~pl.col("barcode").is_in(post["barcode"])
        )
    ).select("barcode", pl.lit(0, dtype=pl.UInt32).alias("n-post-samples"))

    return pl.concat((pre, post))


post_treatment_samples_per_barcode(treated)
# %%
rates = rates.join(post_treatment_samples_per_barcode(treated), on="barcode")

# growth rates vs final culmulative abundance


def percent_vs_rate_group(treated, rates, group):
    rule = alt.Chart(pl.DataFrame({"r": rates["r"].median()})).mark_rule().encode(x="r")
    dots = (
        alt.Chart(treated.filter(pl.col("group") == group).join(rates, on="barcode"))
        .mark_circle(color="black", opacity=0.5)
        .encode(
            alt.X("r")
            .axis(offset=10)
            .title("growth rate (1/hr)")
            .scale(zero=False, domain=[0.015, 0.03]),
            alt.Y("percent").scale(domain=[0, 100]),  # alt.Facet('group').columns(2)
        )
        .properties(title=group, width=200, height=150)
    )
    return rule + dots


@project.save_plot()
def growth_rates_vs_treated_abundance(treated, rates):
    return alt.vconcat(
        alt.hconcat(
            *(percent_vs_rate_group(treated, rates, group) for group in ("PT", "250"))
        ),
        alt.hconcat(
            *(percent_vs_rate_group(treated, rates, group) for group in ("400", "550"))
        ),
    )


growth_rates_vs_treated_abundance(treated, rates)
# %%


# %%
@project.save_plot()
def boxplot_growth_rate_clone_cluster(df, rates):
    width = 40

    return (
        alt.Chart(
            df.join(rates, on="barcode")
            .select("barcode", "clone-cluster", "r")
            .unique()
        )
        .mark_boxplot(
            box=dict(stroke="black", filled=True, size=width),
            rule=dict(
                stroke="black",
            ),
            median=dict(stroke="black", size=width),
            outliers=dict(size=20, filled=True),
            extent="min-max",
        )
        .encode(
            alt.X("clone-cluster").axis(labelAngle=0),
            alt.Y("r").title("growth rate (1/h)").scale(domain=[0.010, 0.03]),
            alt.Color("clone-cluster")
            .scale(**hp.plotting.get_alt_palette("clone-cluster"))
            .legend(None),
        )
        .properties(width=200)
    )


boxplot_growth_rate_clone_cluster(df, rates)

# %%

to_test = df.join(rates, on="barcode").select("barcode", "clone-cluster", "r").unique()
ttest_ind(
    to_test.filter(pl.col("clone-cluster") == "2")["r"],
    to_test.filter(pl.col("clone-cluster") == "1")["r"],
)
# %%
