# %%
import altair as alt
import helpers as hp
import polars as pl
import scanpy as sc
from helpers.plotting.umap import Umap

# %%
project = hp.Project(
    "scrna/downstream", "scrna.downstream.barcodes", outs="outs/new-populations"
)
project.setup_scanpy()
hp.setup_altair()
dataset = hp.data.get_sc_dataset("231-1KB3")
# %%
adata = sc.read(dataset.final)
obs = hp.get_obs(adata)
# %%
# > 20pct of post-sample population
new_populations = [
    ("Dox-250-1", "bg009"),
    ("Dox-250-1", "bg032"),
    ("Dox-250-4", "bg016"),
    ("Dox-250-6", "bg074"),
    ("Dox-250-6", "bg032"),
    ("Dox-400-2", "bg053"),
    ("Dox-400-4", "bg149"),
    ("Dox-400-5", "bg016"),
]


# %%
def load_degs(sample, bgID):
    return pl.read_csv(
        project.paths.data
        / "degs"
        / f"deg-new-population-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    ).with_columns(sample=pl.lit(sample), bgID=pl.lit(bgID))


degs = pl.concat(load_degs(sample, bgID) for sample, bgID in new_populations)


# %%
def get_frequencies(adata, sample, bgID):
    cluster = (
        hp.get_obs(adata)
        .filter(pl.col("sample") == sample, pl.col("bgID") == bgID)
        .select("sample", "bgID", "clone-cluster")
        .unique()["clone-cluster"]
        .item()
    )
    _, var = sc.pp.calculate_qc_metrics(
        adata[(adata.obs["sample"] == sample) & (adata.obs["bgID"] == bgID)],
        inplace=False,
    )

    _, var_ref = sc.pp.calculate_qc_metrics(
        adata[
            (adata.obs["clone-cluster"] == cluster)
            & (adata.obs["sample"] == "Pretreatment")
        ],
        inplace=False,
    )

    return (
        (
            pl.from_pandas(var["pct_dropout_by_counts"].reset_index()).select(
                pl.col("index").alias("gene"),
                (100 - pl.col("pct_dropout_by_counts")).alias("pct_in_group"),
            )
        )
        .join(
            pl.from_pandas(var_ref["pct_dropout_by_counts"].reset_index()).select(
                pl.col("index").alias("gene"),
                (100 - pl.col("pct_dropout_by_counts")).alias("pct_in_ref"),
            ),
            on="gene",
        )
        .with_columns(sample=pl.lit(sample), bgID=pl.lit(bgID))
    )


freqs = pl.concat(
    get_frequencies(adata, sample, bgID) for sample, bgID in new_populations
)
# %%
freqs
# %%
degs = degs.join(
    freqs.filter(pl.col("pct_in_group") > 60, pl.col("pct_in_ref") < 20),
    on=["gene", "sample", "bgID"],
)

# %%
top_genes_each_population = (
    degs.filter(pl.col("FDR") < 0.05, pl.col("pct_in_group") > 0.8)
    .sort(pl.col("log2fc"), descending=True)
    .group_by("bgID", "sample")
    .head()
)
# %%


def top_newly_upregulated(
    degs, fdr_max=0.05, abs_logfc_min=0.25, pct_group=60, pct_reference=25
):
    is_significant_col = (
        (pl.col("FDR") < fdr_max) & (pl.col("log2fc").abs() > abs_logfc_min)
    ).alias("significant")
    return (
        degs.with_columns(
            is_significant_col,
            (pl.col("log2fc").sign())
            .cast(str)
            .replace({"1.0": "up-regulated", "-1.0": "down-regulated"})
            .alias("direction"),
        )
        .filter(
            "significant",
            pl.col("pct_in_group") > pct_group,
            pl.col("pct_in_ref") < pct_reference,
        )
        # .with_columns(
        #     (pl.col("sample") + "-" + pl.col("bgID")).alias("sample-bgID")
        # )
        .filter(pl.col("direction") == "up-regulated")
    )


# %%
newly_upregulated = top_newly_upregulated(degs)


# %%
@project.save_plot()
def bar_newly_upregulated_totals(newly_upregulated):
    source = (
        newly_upregulated.with_columns(
            (pl.col("sample") + "-" + pl.col("bgID")).alias("sample-bgID")
        )
        .group_by("sample-bgID")
        .agg(pl.col("gene").n_unique())
    )
    bar = (
        alt.Chart(source)
        .mark_bar(stroke="#000000", color="#FFFFFF", cornerRadius=2)
        .encode(
            alt.Y("sample-bgID").axis(ticks=False, domain=False, offset=5).title(None),
            alt.X("gene").title("newly upregulated genes (N)"),
        )
        .properties(width=150)
    )
    return bar


bar_newly_upregulated_totals(newly_upregulated)


# %%
def subset_adata(adata_, new_populations):
    labels = (
        hp.get_obs(adata_)
        .filter(pl.col("clone-cluster").is_in(["1", "2"]))
        .with_columns(
            pl.when(pl.col("group").is_in(["250", "400"]))
            .then(pl.col("sample").cast(str) + "-" + pl.col("bgID").cast(str))
            .otherwise(
                pl.col("sample").cast(str)
                + "-clone-cluster-"
                + pl.col("clone-cluster").cast(str)
                + "-barcodes"
            )
            .alias("sample-bgID")
        )
        .select("cell", "sample-bgID")
        .to_pandas()
        .set_index("cell")
    )

    bgIDs = [bgID for _, bgID in new_populations]
    filters = (
        (adata_.obs["n_barcodes"] == 1)
        & (adata_.obs["clone-cluster"].isin(["1", "2"]))
        & (adata_.obs["group"].isin(["PT", "30HR"]) | adata_.obs["bgID"].isin(bgIDs))
    )
    asub = adata_[filters].copy()
    asub.obs = asub.obs.join(labels)
    asub.obs["sample-bgID"] = asub.obs["sample-bgID"].astype("category")
    return asub


asub = subset_adata(adata, new_populations)
# %%

# %%
barcode_newly_upregulated_genes = {
    k[:10]: v
    for k, v in top_newly_upregulated(degs)
    .sort("log2fc", descending=True)
    .filter(~pl.col("gene").str.starts_with("ENSG"))
    .group_by("bgID", maintain_order=True)
    .agg(pl.col("gene").head(10))
    .sort(pl.col("gene").list.len(), pl.col("bgID"), descending=[True, False])
    .rows()
}

# %%
sc.pl.dotplot(
    asub,
    var_names=barcode_newly_upregulated_genes,
    groupby="sample-bgID",
    # dendrogram=True,
    standard_scale="var",
    save="new-upregulated-genes-post-treatment-clones.svg",
)


# %%
def add_sample_bgID_labels(adata_):
    labels = (
        hp.get_obs(adata_)
        .filter(pl.col("clone-cluster").is_in(["1", "2"]))
        .with_columns(
            pl.when(pl.col("group").is_in(["250", "400"]))
            .then(pl.col("sample").cast(str) + "-" + pl.col("bgID").cast(str))
            .otherwise(
                pl.col("sample").cast(str)
                + "-clone-cluster-"
                + pl.col("clone-cluster").cast(str)
                + "-barcodes"
            )
            .alias("sample-bgID")
        )
        .select("cell", "sample-bgID")
        .to_pandas()
        .set_index("cell")
    )

    adata_.obs = adata_.obs.join(labels.astype("category"))
    # adata_.obs['sample-bgID'] = adata_.obs['sample-bgID'].astype('category')


# %%
add_sample_bgID_labels(adata)


# %%
def new_populations_umap(
    adata_,
    sample,
    bgID,
):
    adata_ = adata_.copy()
    obs = hp.get_obs(adata_)
    cluster = obs.filter(pl.col("bgID") == bgID)["clone-cluster"].unique().item()

    cells = (
        hp.get_obs(adata_)
        .filter(
            pl.col("sample-bgID").is_in(
                [f"Pretreatment-clone-cluster-{cluster}-barcodes", f"{sample}-{bgID}"]
            )
        )["cell"]
        .to_list()
    )
    Umap(adata_, "sample-bgID", subset_idx=cells).plot(
        cmap={
            f"Pretreatment-clone-cluster-{cluster}-barcodes": hp.plotting.get_palette(
                "sample"
            )["Pretreatment"],
            f"{sample}-{bgID}": hp.plotting.get_palette("sample")[sample],
        },
        legend_kwargs=dict(
            bbox_to_anchor=(0.04, 1),
            loc="center left",
        ),
        title="",
        save=project.paths.outs / f"umap-{sample}-{bgID}-Pretreatment-cluster.png",
    )


# %%

for sample, bgID in new_populations:
    new_populations_umap(adata, sample, bgID)
# %%
# ---
