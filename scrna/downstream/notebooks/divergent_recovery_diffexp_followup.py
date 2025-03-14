# %%
import sys
from textwrap import wrap

import altair as alt
import helpers as hp
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from helpers.scrna.gsea import split_gsea_term_column

# %%
project = hp.Project(
    "scrna/downstream",
    "scrna.downstream.barcodes",
    outs="outs/divergent",
)
project.setup_scanpy()
hp.setup_altair()

divergent_bgIDs_samples = [
    ("bg032", "Dox-250-1"),
    ("bg032", "Dox-250-6"),
    ("bg016", "Dox-250-4"),
    ("bg016", "Dox-400-5"),
    ("bg010", "Dox-250-4"),
    ("bg010", "Dox-400-4"),
]
divergent_barcodes = list(set(bgID for bgID, _ in divergent_bgIDs_samples))


# %%
def load_degs(bgID, sample):
    return pl.read_csv(
        project.paths.data
        / "degs"
        / f"deg-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    ).with_columns(bgID=pl.lit(bgID), sample=pl.lit(sample))


degs = pl.concat((load_degs(bgID, sample) for bgID, sample in divergent_bgIDs_samples))


# %%
def add_significance(degs):
    return degs.with_columns(
        pl.when(pl.col("log2fc").sign() > 0)
        .then(pl.lit("up"))
        .otherwise(pl.lit("down"))
        .alias("direction"),
        ((pl.col("FDR") < 0.1) & (pl.col("log2fc").abs() > 0.25)).alias("significant"),
    )


degs = add_significance(degs)
# %%
# Total differentially expressed genes


def overlapping_genes(degs, bgID):
    source = degs.filter((pl.col("bgID") == bgID) & (pl.col("significant")))
    source = (
        source.join(
            source.group_by("gene", "direction").agg(
                pl.col("sample").n_unique().alias("n_samples")
            ),
            on="gene",
        ).with_columns(
            pl.when(pl.col("n_samples") > 1)
            .then(pl.lit("both"))
            .otherwise(pl.col("sample"))
            .alias("sample")
        )
        # .select("gene", "direction", "sample")
        # .unique()
    )
    source = (
        source.group_by("sample", "direction")
        .len()
        .with_columns(pl.col("len").alias("n_genes"))
    )
    charts = []
    max_x = source["n_genes"].max() * 1.1
    sample_list = [
        "both",
        *sorted(sample for sample in source["sample"].unique() if sample != "both"),
    ]
    source = source.with_columns(pl.col("direction") + "-regulated")

    colors = {"blue": "#0000FF", "grey": "#808080", "red": "#FF0000"}
    color_scale = dict(
        domain=["up-regulated", "down-regulated"], range=[colors["blue"], colors["red"]]
    )

    # TODO: use different values for color scale
    for i, sample in enumerate(sample_list):
        x = alt.X("n_genes").title("sig. DEGs (N)").scale(domain=[0, max_x])
        if i != 2:
            x = x.axis(None)

        base = (
            alt.Chart(source.filter(pl.col("sample") == sample))
            .mark_bar()
            .encode(
                x,
                alt.Y("direction")
                # hack to add sample labels
                .title(sample)
                .axis(
                    domain=False,
                    ticks=False,
                    labels=False,
                    titleAngle=0,
                    titleAlign="right",
                    titleBaseline="middle",
                    titlePadding=20,
                ),
            )
            .properties(width=150, height=60)
        )

        bar = base.encode(
            alt.Color("direction").scale(**color_scale).title(None),
        )
        text = base.mark_text(color="black", dx=15).encode(alt.Text("n_genes"))
        charts.append(bar + text)
    return alt.vconcat(*charts, spacing=2).configure_legend(
        orient="none", legendX=100, legendY=0
    )


# %%
for bgID in divergent_barcodes:
    p = overlapping_genes(degs, bgID)
    p.save(project.paths.outs / f"overlapping-genes-{bgID}-pre-vs-post-treatment.svg")
    p.show()


# %%
def volcano_plot(dedf, pval_cutoff=0.05, logFC_cutoff=0.25):
    if dedf["pval"].min() == 0:
        dedf = dedf.with_columns(pl.col("pval") + sys.float_info.min)
    fig, ax = plt.subplots()

    line_kwargs = dict(color="black", linestyle="--")

    sig_genes = pl.col("pval") < pval_cutoff
    up_genes = (pl.col("log2fc") > logFC_cutoff) & sig_genes
    down_genes = (pl.col("log2fc") < -1 * logFC_cutoff) & sig_genes
    else_genes = ~(up_genes | down_genes)

    for f, color in ((up_genes, "blue"), (down_genes, "red"), (else_genes, "grey")):
        ax.scatter(
            x=dedf.filter(f)["log2fc"],
            y=dedf.filter(f)["pval"].log10() * -1,
            c=color,
            s=5,
            alpha=0.8,
        )
    ax.annotate(
        f"{dedf.filter(up_genes).shape[0]} up-regulated",
        (0.95, 1.01),
        xycoords="axes fraction",
        ha="right",
    )
    ax.annotate(
        f"{dedf.filter(down_genes).shape[0]} down-regulated",
        (0.05, 1.01),
        xycoords="axes fraction",
        ha="left",
    )

    x_lim = dedf["log2fc"].abs().max() * 1.1
    ax.set_xlim(-1 * x_lim, x_lim)
    ax.axhline(-1 * np.log10(pval_cutoff), **line_kwargs)
    ax.axvline(logFC_cutoff, **line_kwargs)
    ax.axvline(-1 * logFC_cutoff, **line_kwargs)
    ax.set(xlabel=r"$\log_{2}(FC)$", ylabel=r"$-\log_{10}(pval)$")
    return fig


# %%
def plot_both_volcanos(degs, bgID):
    degs = degs.filter(pl.col("bgID") == bgID)
    samples = degs["sample"].unique().to_list()
    for sample in samples:
        volcano_plot(degs.filter(pl.col("sample") == sample))
        plt.savefig(
            project.paths.outs / f"divergent-barcode-{barcode}-{sample}-deg-volcano.png"
        )


for barcode in divergent_barcodes:
    plot_both_volcanos(degs, barcode)


# %%
def load_gsea(bgID, sample):
    return pl.read_csv(
        project.paths.data
        / "gsea"
        / f"gsea-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    ).with_columns(pl.lit(sample).alias("sample"), pl.lit(bgID).alias("bgID"))


gsea = pl.concat((load_gsea(bgID, sample) for bgID, sample in divergent_bgIDs_samples))
gsea = split_gsea_term_column(gsea)


# %%
def plot_gsea(gsea):
    source = gsea.with_columns(
        (-1 * (pl.col("FDR q-val").replace([0], 1 / 10000)).log10()).alias(
            "-log10(FDR)"
        )
    )
    color_scale = hp.plotting.get_alt_palette(
        "sample", source["sample"].unique().to_list()
    )
    source = source.with_columns(
        pl.col("term")
        .replace_strict(
            {term: "|".join(wrap(term, width=35)) for term in gsea["term"].unique()}
        )
        .alias("y")
    )
    term_order = (
        source.group_by("y").agg(pl.col("NES").max()).sort("NES", descending=True)["y"]
    )
    return (
        alt.Chart(source)
        .mark_circle()
        .encode(
            alt.X("NES"),
            alt.Y("y")
            .title("GO Biological Process Terms")
            .scale(domain=term_order)
            .axis(offset=10),
            alt.Size("-log10(FDR)"),
            alt.Color("sample").scale(**color_scale),
        )
        # has to be done at the label level or it breaks the chart
        # can't find any existing altair bug report besides this stack overflow:
        # labelExpr: https://stackoverflow.com/a/77456810
        # labelOffset: https://stackoverflow.com/a/77548746
        .configure_axisX(grid=True)
        .configure_axisY(
            labelExpr='split(datum.label, "|")',
            labelLineHeight=10,
            labelOffset=alt.ExprRef(
                'length(split(datum.label, "|")) > 1 ? (length(split(datum.label, "|")) > 2 ? -10 : -5) : 0'
            ),
            labelLimit=200,
        )
        .properties(height=250, width=100)
    )


def plot_gsea_barcode(gsea, bgID):
    gsea = gsea.filter(pl.col("bgID") == bgID)
    top = (
        gsea.filter(pl.col("FDR q-val") < 0.05)
        .sort("FDR q-val", pl.col("NES").abs(), descending=[False, True])
        .group_by("sample", maintain_order=True)
        .head(5)
    )
    p = plot_gsea(gsea.filter(pl.col("term").is_in(top["term"])))
    p.show()
    p.save(project.paths.outs / f"gsea-top5-significant-{bgID}.svg")


# %%

for bgID in divergent_barcodes:
    plot_gsea_barcode(
        gsea,
        bgID,
    )
# %%
# ---
