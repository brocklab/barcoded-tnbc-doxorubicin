# %%
import helpers as hp
import polars as pl
import scanpy as sc

# %%
project = hp.Project(
    "scrna/downstream", "scrna.downstream", outs="outs/deg-clone-cluster"
)
project.setup_scanpy()
hp.setup_altair()
dataset = hp.data.get_sc_dataset("231-1KB3")
# %%
adata = sc.read(dataset.final)


# %%
def load_degs(sample):
    return pl.read_csv(
        project.paths.data / "degs" / f"deg-{sample}-clone-cluster-1-vs-2.csv"
    ).with_columns(sample=pl.lit(sample))


degs = pl.concat(load_degs(sample) for sample in ["Pretreatment", "Dox-550-30HR"])

# %%
adata_sub = adata[adata.obs["clone-cluster"].isin(["1", "2"])].copy()
# %%
adata_sub.obs["sample/clone-cluster"] = (
    adata_sub.obs["sample"].astype("str")
    + "/"
    + adata_sub.obs["clone-cluster"].astype("str")
).astype("category")

# one cell
adata_sub = adata_sub[adata_sub.obs["sample/clone-cluster"] != "Dox-400-2/2"].copy()


adata_sub.obs["sample/clone-cluster"] = (
    adata_sub.obs["clone-cluster"].map({"1": "", "2": " "}).astype(str)
    + adata.obs["sample"].astype(str)
).astype("category")

adata_sub.obs["sample/clone-cluster"] = adata_sub.obs[
    "sample/clone-cluster"
].cat.set_categories(
    list(
        dict.fromkeys(
            adata_sub.obs.sort_values(by=["clone-cluster", "sample"])[
                "sample/clone-cluster"
            ]
        )
    )
)

adata_sub.uns["sample/clone-cluster_colors"] = [
    hp.plotting.get_palette("sample")[k.strip()]
    for k in adata_sub.obs["sample/clone-cluster"].cat.categories
]


# %%
def genes_for_heatmap(degs, n_genes=500):
    genes = (
        degs.filter(
            pl.col("log2fc").abs() > 0.25,
            pl.col("FDR") < 0.05,
            pl.col("sample") == "Pretreatment",
        )
        .sort("log2fc", descending=True)["gene"]
        .to_list()
    )
    assert len(genes) > n_genes
    return {
        "clone-cluster 1 markers": genes[-int(n_genes / 2) :],
        "clone-cluster 2 markers": genes[: int(n_genes / 2)],
    }


# %%
hp.plotting.heatmap(
    adata_sub,
    genes=genes_for_heatmap(degs, n_genes=400),
    groupby="sample/clone-cluster",
    var_group_rotation=0,
    figsize=(8, 5),
    save=project.paths.outs / "heatmap-pretreament-clone-cluster-1-vs-2.svg",
)

# %%
