# TODO: cleanup
# %%
from collections import namedtuple

import altair as alt
import helpers as hp
import numpy as np
import polars as pl
from scipy.optimize import curve_fit

# %%
project = hp.Project("drug-dosing")
hp.setup_altair()
# %%
# plate -> well -> sample
platemaps = dict(
    DM2470A={
        **{f"{c}{num}": "MDA-MB-231" for num in range(1, 13) for c in "ABCD"},
        **{f"{c}{num}": "231-1KB3" for num in range(1, 13) for c in "EFGH"},
        **{"A12": "blank", "E12": "blank"},
    },
    DM2470B={
        **{f"{c}{num}": "231-1KB3-EN2" for num in range(1, 13) for c in "ABCD"},
        **{f"{c}{num}": "231-1KB3-EP2" for num in range(1, 13) for c in "EFGH"},
        **{"A12": "blank", "E12": "blank"},
    },
)
samples = ["MDA-MB-231", "231-1KB3", "231-1KB3-EN2", "231-1KB3-EP2"]

# column -> dose
doses = {
    # max dose 10uM decreasing by a factor of 3 starting from column 2
    **{i + 2: 10000 / (3**i) for i in range(9)},
    # control
    **{11: 0},
}
inner60_wells = [f"{c}{num}" for c in "BCDEFG" for num in range(2, 12)]


# %%
# def read_confluence_data(plate, platemap):
#     assert plate in platemap, "expected plate to exist in platemap"
#     f = project.paths.data / f"{plate}-incucyte-confluence.txt"
#     assert f.is_file(), f"expected {f} to exist"
#     df = pl.read_csv(f, separator="\t", skip_rows=1).rename(
#         {"Date Time": "time", "Elapsed": "elapsed"}
#     )

#     return (
#         df.unpivot(
#             index=["time", "elapsed"], variable_name="well", value_name="percent"
#         )
#         .with_columns(
#             # one of them is int, so let's cast then both to float
#             pl.col("elapsed").cast(float),
#             pl.col("well").replace_strict(platemaps[plate]).alias("sample"),
#             pl.col("well").str.slice(0, 1).alias("row"),
#             pl.col("well").str.slice(1).cast(int).alias("column"),
#             pl.lit(plate).alias("plate"),
#         )
#         .with_columns()
#     )


# df = pl.concat(read_confluence_data(p, platemaps) for p in platemaps).filter(
#     pl.col("sample") != "blank"
# )


# # %%
# def microplate(df, plate):
#     base = (
#         alt.Chart(
#             df.filter(pl.col("plate") == plate, pl.col("well").is_in(inner60_wells)),
#         )
#         .encode(
#             alt.X("elapsed").title(None).scale(padding=3).axis(domain=False),
#             alt.Y("percent").title(None).axis(offset=-1, domain=False).scale(padding=3),
#         )
#         .properties(width=50, height=50)
#     )
#     return (
#         (base.mark_point() + base.mark_line())
#         .facet(
#             row=alt.Row("row").title("percent").header(labelAngle=0),
#             column=alt.Column("column:O")
#             .title("elapsed time")
#             .header(titleOrient="bottom"),
#         )
#         .properties(spacing=0, title=plate)
#         # .configure_view(strokeWidth=1,stroke='black')
#     )


# (
#     alt.vconcat(microplate(df, "DM2470A"), microplate(df, "DM2470B"))
#     .configure_view(strokeWidth=1, stroke="black")
#     .configure_title(align="center", anchor="middle")
# )
# # %%


# def sample_by_dose(df, sample):
#     source = (
#         df.filter(pl.col("sample") == sample, pl.col("well").is_in(inner60_wells))
#         .with_columns(pl.col("column").replace_strict(doses).alias("dose"))
#         .group_by("sample", "dose", "elapsed")
#         .agg(
#             pl.col("percent").mean(),
#             (pl.col("percent").std() / math.sqrt(3)).alias("sem"),
#         )
#     ).with_columns(
#         (pl.col("percent") - pl.col("sem")).alias("sem-low"),
#         (pl.col("percent") + pl.col("sem")).alias("sem-high"),
#     )
#     base = alt.Chart(source).encode(
#         alt.X("elapsed").title("elapsed time"), alt.Y("percent"), alt.Color("dose:O")
#     )

#     errorbars = base.encode(alt.Y("sem-low"), alt.Y2("sem-high")).mark_errorbar(
#         ticks=True
#     )
#     return (base.mark_circle() + base.mark_line() + errorbars).properties(title=sample)


# (
#     (sample_by_dose(df, "MDA-MB-231") | sample_by_dose(df, "231-1KB3"))
#     & (sample_by_dose(df, "231-1KB3-EN2") | sample_by_dose(df, "231-1KB3-EP2"))
# )
# # %%
# doses


# %%
# determine change from "initial time", loop per sample

# def percent_change(df, sample):
#     cols = ["sample", "dose"]
#     df = df.filter(
#         pl.col("sample") == sample, pl.col("well").is_in(inner60_wells)
#     ).with_columns(pl.col("column").replace_strict(doses).alias("dose"))

#     return (
#         (
#             df.filter(pl.col("elapsed") == df["elapsed"].min())
#             .select((*cols, "percent"))
#             .rename({"percent": "initial"})
#         )
#         .join(
#             df.filter(pl.col("elapsed") == df["elapsed"].max())
#             .select((*cols, "percent"))
#             .rename({"percent": "final"}),
#             on=cols,
#         )
#         .with_columns(
#             ((pl.col("final") - pl.col("initial")) / pl.col("final") * 100).alias(
#                 "percent change"
#             )
#         )
#     )


# def percent_change_bar(df):
#     source = (
#         pl.concat(percent_change(df, sample) for sample in df["sample"].unique())
#         .group_by("sample", "dose")
#         .agg(pl.col("percent change").mean())
#     )

#     return (
#         alt.Chart(source)
#         .mark_bar()
#         .encode(
#             alt.Y("dose:O"), alt.X("percent change"), alt.Facet("sample").columns(2)
#         )
#     )


# percent_change_bar(df)


# %%
def read_celltiter_data(plate, platemap):
    assert plate in platemap, "expected plate to exist in platemap"
    f = project.paths.data / f"{plate}-celltiter.txt"
    assert f.is_file(), f"expected {f} to exist"
    return (
        pl.read_csv(
            project.paths.data / f"{plate}-celltiter.txt",
            skip_rows=1,
            has_header=False,
            new_columns=["row", *map(str, range(1, 13))],
        )
        .unpivot(index=["row"], variable_name="column", value_name="lum-raw")
        .filter(pl.col("lum-raw").is_not_null())
        .with_columns(
            pl.col("lum-raw").cast(int),
            (pl.col("row") + pl.col("column")).alias("well"),
            pl.col("column").cast(int),
        )
        .with_columns(
            pl.col("well").replace_strict(platemaps[plate]).alias("sample"),
            pl.lit(plate).alias("plate"),
        )
    )


def process_dosing_data(df):
    blanks = dict(
        df.filter(pl.col("sample") == "blank")
        .group_by("plate")
        .agg(pl.col("lum-raw").mean())
        .rows()
    )

    df = (
        df.filter(pl.col("sample") != "blank")
        .with_columns(
            (pl.col("lum-raw") - pl.col("plate").replace_strict(blanks)).alias("lum")
        )
        .filter(pl.col("well").is_in(inner60_wells))
        .with_columns(pl.col("column").replace_strict(doses).alias("dose"))
    )

    control_lum = dict(
        df.filter(pl.col("dose") == 0)
        .group_by("sample")
        .agg(pl.col("lum").mean())
        .rows()
    )
    df = df.filter(pl.col("dose") != 0).with_columns(
        (pl.col("lum") / pl.col("sample").replace_strict(control_lum) * 100).alias(
            "viability"
        )
    )
    return df


df = process_dosing_data(
    pl.concat(read_celltiter_data(p, platemaps) for p in platemaps)
)


# %%
class Fit:
    def __init__(self, df, func, params_nt, p0, bounds=(-np.inf, np.inf)):
        self.df = df
        self.func = func
        self.params_nt = params_nt
        self.p0 = p0
        self.bounds = bounds
        self._apply_fit()

    def _apply_fit(self):
        x = self.df["dose"]
        y = self.df["viability"]
        params, cov = curve_fit(self.func, x, y, p0=self.p0, bounds=self.bounds)
        self.params = self.params_nt(*params)
        self.cov = cov

    def __repr__(self):
        return f"Fit(func={self.func.__name__}, {self.params})"

    def get_df(self):
        x_fit = np.logspace(0, 5, 500)
        return pl.DataFrame(
            {"dose": x_fit, "viability": self.func(x_fit, *self.params)}
        )


def fourPL(x, bottom, top, hill, ec50):
    """four-parameter logistic equation a.k.a. hill equation
    with constrained top and bottom
    """
    return bottom + (top - bottom) / (1 + (x / ec50) ** hill)


unconstrained_fourPL_fits = {
    sample: Fit(
        df.filter(pl.col("sample") == sample),
        fourPL,
        namedtuple("params", "bottom top hill IC50"),
        p0=[0.0, 100, 1.0, 300],
    )
    for sample in samples
}


def constrained_fourPL(x, hill, ec50):
    """four-parameter logistic equation a.k.a. hill equation
    with constrained top and bottom
    """
    bottom = 0
    top = 100
    return bottom + (top - bottom) / (1 + (x / ec50) ** hill)


fourPL_fits = {
    sample: Fit(
        df.filter(pl.col("sample") == sample),
        constrained_fourPL,
        namedtuple("params", "hill IC50"),
        p0=[1.0, 300],
    )
    for sample in samples
}


# %%
def drug_dose_fit_curve(sample, fit):
    pts = (
        alt.Chart(fit.df)
        .mark_point()
        .encode(
            alt.X("dose")
            .title("Doxorubicin (nM)")
            .scale(type="log")
            .axis(
                # hack to make it show only the power of 10 labels
                labelExpr="datum.label[0] == '1' ? datum.label : ''"
            ),
            alt.Y("viability").title("percent viability"),
        )
    )
    line = alt.Chart(fit.get_df()).mark_line().encode(alt.X("dose"), alt.Y("viability"))
    label = (
        alt.Chart()
        .mark_text(x=150, y=25)
        .encode(text=alt.value(f"IC50 = {fit.params.IC50:.2f} nM"))
    )
    return (pts + line + label).properties(title=sample, width=200)


def plot_combined_dose_response_curves(fits):
    panels = {sample: drug_dose_fit_curve(sample, fits[sample]) for sample in samples}
    return (panels["MDA-MB-231"] | panels["231-1KB3"]) & (
        panels["231-1KB3-EN2"] | panels["231-1KB3-EP2"]
    )


# %%
@project.save_plot()
def combined_dose_response_curves_unconstrained():
    return plot_combined_dose_response_curves(unconstrained_fourPL_fits)


combined_dose_response_curves_unconstrained()

# %%

for sample, fit in unconstrained_fourPL_fits.items():
    p = drug_dose_fit_curve(sample, fit)
    p.save(project.paths.outs / f"dose-response-curve-{sample}.svg")
    p.show()


# %%
@project.save_plot()
def combined_dose_response_curves():
    return plot_combined_dose_response_curves(fourPL_fits)


combined_dose_response_curves()


# %%
def single_curve(sample, fit):
    return (
        alt.Chart(fit.get_df().with_columns(sample=pl.lit(sample)))
        .mark_line()
        .encode(
            alt.X("dose")
            .title("Doxorubicin (nM)")
            .scale(type="log")
            .axis(
                # hack to make it show only the power of 10 labels
                labelExpr="datum.label[0] == '1' ? datum.label : ''"
            ),
            alt.Y("viability").title("percent viability"),
            alt.Color("sample").sort(samples),
        )
    )


def combined_curves(fits):
    curves = [single_curve(sample, fit) for sample, fit in fits.items()]
    return alt.layer(*curves).configure_axis(offset=5)


# %%
@project.save_plot()
def dox_response_all_samples_curves_only_unconstrained():
    return combined_curves(unconstrained_fourPL_fits)


dox_response_all_samples_curves_only_unconstrained()


# %%
@project.save_plot()
def dox_response_all_samples_curves_only():
    return combined_curves(fourPL_fits)


dox_response_all_samples_curves_only()


# %%
@project.save_plot()
def IC50_bar(fits):
    source = pl.DataFrame(
        ((sample, fit.params.IC50) for sample, fit in fits.items()),
        schema=("sample", "IC50"),
    )
    bar = (
        alt.Chart(source)
        .mark_bar(cornerRadius=3, stroke="black", color="white")
        .encode(
            alt.Y("sample").sort(samples).axis(domain=False, ticks=False, offset=2),
            alt.X("IC50").title("Doxorubicin IC50 (nM)").axis(offset=2),
        )
    )
    text = bar.mark_text(xOffset=25).encode(alt.Text("IC50").format(".2f"))
    return bar + text


IC50_bar(fourPL_fits)
# %%
