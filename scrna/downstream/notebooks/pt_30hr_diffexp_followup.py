# %%
from textwrap import wrap

import altair as alt
import helpers as hp
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scanpy as sc
from helpers.plotting import volcano
from helpers.scrna.gsea import split_gsea_term_column
from scipy.stats import linregress

alt.data_transformers.enable("vegafusion")

# %%
project = hp.Project("scrna/downstream", "scrna.downstream", outs="outs/pt-30hr")
project.setup_scanpy()
hp.setup_altair()
dataset = hp.data.get_sc_dataset("231-1KB3")
# %%
adata = sc.read(dataset.final)
# %%
# DEG analysis only considered clones with known identity
adata = adata[adata.obs["clone-cluster"].isin(["1", "2"])].copy()
adata.obs["sample/clone-cluster"] = (
    adata.obs["sample"].astype("str") + "/" + adata.obs["clone-cluster"].astype("str")
).astype("category")


def load_degs(cluster):
    return (
        pl.read_csv(
            project.paths.data
            / "degs"
            / f"deg-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
        )
    ).with_columns(cluster=pl.lit(cluster))


degs = pl.concat((load_degs(cluster) for cluster in ("1", "2")))
degs.head()
# %%
print("number of significant genes found in both groups:")
degs.filter(pl.col("FDR") < 0.05).group_by("gene").len()["len"].value_counts()


# %%
def pivot_degs(degs):
    overlap = degs.filter(
        pl.col("gene").is_in(
            degs.group_by("gene").len().filter(pl.col("len") == 2)["gene"]
        )
    )
    return overlap.pivot(on="cluster", index="gene", values=["log2fc"])


def fit_line(df):
    slope, intercept, rvalue, _, _ = linregress(df["1"], df["2"])
    return rvalue**2, pl.DataFrame(
        {"x": (x := np.linspace(-2.5, 2.5)), "y": (slope * x + intercept)}
    )


@project.save_plot()
def pt_vs_30hr_clone_cluster_deg_correlation(degs):
    source = pivot_degs(degs)
    rsquared, line_df = fit_line(source)

    line = (
        alt.Chart(line_df)
        .mark_line(
            color="red",
        )
        .encode(alt.X("x"), alt.Y("y"))
    )

    def gen_enc(encoding, cluster):
        return encoding.scale(domain=[-3, 3]).title(
            f"log₂(fold change) clone-cluster {cluster}"
        )

    label = alt.Chart().mark_text(
        x=10, y=5, text=f"R²={rsquared:.2f}", align="left", fontStyle="bold"
    )
    points = (
        alt.Chart(source)
        .mark_circle(size=3, opacity=0.5, color="black")
        .encode(gen_enc(alt.X("1"), "1"), gen_enc(alt.Y("2"), "2"))
    )
    return points + line + label


pt_vs_30hr_clone_cluster_deg_correlation(degs)


# %%
def load_gsea(cluster):
    return pl.read_csv(
        project.paths.data
        / "gsea"
        / f"gsea-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
    ).with_columns(cluster=pl.lit(cluster))


gsea = pl.concat([load_gsea(cluster) for cluster in ("1", "2")])
gsea = split_gsea_term_column(gsea)

# %%
# sig_one_terms = gsea.filter(pl.col("FDR q-val") < 0.05)["term"].unique().to_list()
# print(sig_one_terms)
# list paired down from sig_one_terms
curated_terms = [
    "Positive Regulation Of Cellular Metabolic Process",
    "tRNA Aminoacylation For Protein Translation",
    "Mismatch Repair",
    "Antigen Processing And Presentation Of Peptide Antigen Via MHC Class Ib",
    "DNA Replication",
    "ATP Biosynthetic Process",
    "Negative Regulation Of mRNA Splicing, Via Spliceosome",
    "DNA Metabolic Process",
    "NADH Dehydrogenase Complex Assembly",
    "DNA Recombination",
    "Cotranslational Protein Targeting To Membrane",
    "Cytoplasmic Translation",
    "DNA Replication Initiation",
    "'De Novo' Post-Translational Protein Folding",
]
# %%


@project.save_plot()
def pt_vs_30hr_clone_clusters_gsea(gsea):
    color_scale = hp.plotting.get_alt_palette("clone-cluster", ("1", "2"))

    source = gsea.filter(pl.col("term").is_in(curated_terms)).with_columns(
        (pl.col("FDR q-val").log10() * -1).alias("-log10(FDR)")
    )

    source = source.with_columns(
        pl.col("term")
        .replace_strict(
            {term: "|".join(wrap(term, width=35)) for term in gsea["term"].unique()}
        )
        .alias("y")
    )
    term_order = (
        source.group_by("y")
        .agg(pl.col("NES").max())
        .sort("NES", descending=True)["y"]
        .to_list()
    )

    print(term_order)

    return (
        alt.Chart(source)
        .mark_circle()
        .encode(
            alt.X("NES"),
            alt.Y("y").title("GO Biological Process Terms").scale(domain=term_order),
            alt.Size("-log10(FDR)"),
            alt.Color("cluster").scale(**color_scale),
        )
        # has to be done at the label level or it breaks the chart
        # can't find any existing altair bug report besides this stack overflow:
        # labelExpr: https://stackoverflow.com/a/77456810
        # labelOffset: https://stackoverflow.com/a/77548746
        .configure_axisX(grid=True)
        .configure_axisY(
            labelExpr='split(datum.label, "|")',
            labelLineHeight=10,
            labelOffset=alt.ExprRef('length(split(datum.label, "|")) != 1 ? -5 : 0'),
            labelLimit=200,
        )
        .properties(height=325, width=150)
    )


pt_vs_30hr_clone_clusters_gsea(gsea)


# %%
def plot_heatmap(adata, degs, gene_cutoff=500, **kwargs):
    samples = (
        hp.get_obs(adata)
        .group_by("sample/clone-cluster")
        .len()
        .filter(pl.col("len") > 20)
    )["sample/clone-cluster"].to_list()

    # filter by samples with enough cells
    adata = adata[adata.obs["sample/clone-cluster"].isin(samples)].copy()

    genes_by_log2fc = (
        degs.group_by("gene")
        .agg(pl.col("log2fc").mean())
        .sort("log2fc", descending=True)["gene"]
        .to_list()
    )
    genes = {
        "up-regulated": genes_by_log2fc[:gene_cutoff],
        "down-regulated": genes_by_log2fc[-gene_cutoff:],
    }
    # all_genes = (genes_by_log2fc[:gene_cutoff] + genes_by_log2fc[-gene_cutoff:],)

    # recalculate dendrogram with important genes?

    # sc.tl.dendrogram(adata, groupby="sample/clone-cluster", var_names=all_genes)

    hp.plotting.heatmap(
        adata,
        genes=genes,
        groupby="sample/clone-cluster",
        var_group_rotation=0,
        **kwargs,
    )


# %%
adata_hmap = adata[adata.obs["clone-cluster"].isin(["1", "2"])].copy()
adata_hmap.obs["sample/clone-cluster"] = (
    adata_hmap.obs["sample"].astype("str")
    + "/"
    + adata_hmap.obs["clone-cluster"].astype("str")
).astype("category")

# one cell
adata_hmap = adata_hmap[adata_hmap.obs["sample/clone-cluster"] != "Dox-400-2/2"].copy()


adata_hmap.obs["sample/clone-cluster"] = (
    adata_hmap.obs["clone-cluster"].map({"1": "", "2": " "}).astype(str)
    + adata.obs["sample"].astype(str)
).astype("category")

adata_hmap.obs["sample/clone-cluster"] = adata_hmap.obs[
    "sample/clone-cluster"
].cat.set_categories(
    list(
        dict.fromkeys(
            adata_hmap.obs.sort_values(by=["clone-cluster", "sample"])[
                "sample/clone-cluster"
            ]
        )
    )
)

adata_hmap.uns["sample/clone-cluster_colors"] = [
    hp.plotting.get_palette("sample")[k.strip()]
    for k in adata_hmap.obs["sample/clone-cluster"].cat.categories
]


# %%

plot_heatmap(
    adata_hmap,
    degs,
    gene_cutoff=250,
    save=project.paths.outs / "heatmap-pt-vs-30hr-clone-clusters-deg-up250-down250.svg",
)

# %%
for cluster in "1", "2":
    volcano(
        degs.filter(pl.col("cluster") == cluster),
        save=project.paths.outs / f"pt-vs-30hr-clone-cluster-{cluster}-deg-volcano.png",
        save_kwargs=dict(dpi=300),
    )
    plt.show()

# %%
