# TODO: use config.yml for all colors, fix survivorship
# TODO: no cell with barcode should have "unknown" survivorship
# %%
import helpers as hp
import polars as pl
import scanpy as sc

# %%
project = hp.Project("scrna/downstream")
project.setup_scanpy()
# %%
dataset = hp.data.get_sc_dataset("231-1KB3")
adata = sc.read_h5ad(dataset.barcoded)
targeted = hp.data.targeted_final.load()

# %%
# set default category and order
adata.obs["group"] = adata.obs["group"].cat.set_categories(["PT", "30HR", "250", "400"])
adata.obs["sample"] = adata.obs["sample"].cat.set_categories(
    [
        "Pretreatment",
        "Dox-550-30HR",
        "Dox-250-1",
        "Dox-250-4",
        "Dox-250-6",
        "Dox-400-2",
        "Dox-400-4",
        "Dox-400-5",
    ]
)


# %%
def survivorship_from_targeted(targeted):
    bar_known = pl.col("barcode") != "unknown"
    is_pt = pl.col("sample") == "PT-0"
    pre_barcodes = targeted.filter(bar_known, is_pt).select("barcode")
    post_barcodes = targeted.filter(bar_known, ~is_pt).select("barcode")
    in_pre_barcodes = pl.col("barcode").is_in(pre_barcodes)
    in_post_barcodes = pl.col("barcode").is_in(post_barcodes)

    post_targeted_barcodes_survivorship = (
        targeted.filter(in_pre_barcodes, in_post_barcodes)
        .join(
            targeted.filter(bar_known, is_pt).select("barcode", pt_percent="percent"),
            on="barcode",
        )
        .with_columns(
            (
                ((pl.col("percent") - pl.col("pt_percent")) / pl.col("pt_percent"))
                * 100
            ).alias("change")
        )
        .group_by("barcode")
        .agg(pl.col("change").max())
        .with_columns(
            pl.when(pl.col("change") > 5)
            .then(pl.lit("high"))
            .otherwise(pl.lit("low"))
            .alias("survivorship")
        )
    ).select("barcode", "survivorship")

    return pl.concat(
        [
            post_targeted_barcodes_survivorship,
            targeted.filter(in_pre_barcodes, ~in_post_barcodes)
            .with_columns(pl.lit("low").alias("survivorship"))
            .select("barcode", "survivorship"),
            targeted.filter(in_post_barcodes, ~in_pre_barcodes)
            .with_columns(pl.lit("high").alias("survivorship"))
            .select("barcode", "survivorship")
            .unique(),
        ]
    )


adata.obs["survivorship"] = (
    adata.obs["barcode"]
    .map(
        survivorship_from_targeted(targeted)
        .to_pandas()
        .set_index("barcode")["survivorship"]
    )
    .fillna("unknown")
    .astype("category")
)

# barcodes only usually found in the library are also low survivorship
adata.obs.loc[
    (adata.obs["survivorship"] == "unknown") & (adata.obs["n_barcodes"] == 1),
    "survivorship",
] = "low"

# %%
bgIDs = dict(hp.data.barcode_bgID.load().rows())


# generate a mapping for all barcodes + barcode-pairs
def barcodes_to_bgIDs(barcodes):
    return "-".join((bgIDs[bar] for bar in barcodes))


def map_barcodes_to_bgIDs(adata):
    return dict(
        hp.get_obs(adata)
        .filter(pl.col("barcode").is_not_null())
        .select("barcode")
        .with_columns(
            pl.col("barcode")
            .cast(str)
            .str.split("-")
            .map_elements(barcodes_to_bgIDs, return_dtype=str)
            .alias("bgID")
        )
        .rows()
    )


adata.obs["bgID"] = adata.obs["barcode"].map(map_barcodes_to_bgIDs(adata))

# %%
adata.obs["clone-specific"] = adata.obs["n_barcodes"] == 1


# %%
def add_palette_to_uns(adata, category, categories=None):
    palette = hp.plotting.get_palette(category, categories)
    items = adata.obs[category].cat.categories.to_list()
    assert set(items) == set(palette.keys()), (
        f"\ncategories do not match\n\tobs: {sorted(items)}\n\tpalette: {sorted(palette.keys())} "
    )
    adata.uns[f"{category}_colors"] = [palette[x] for x in items]


# %%
add_palette_to_uns(adata, "clone-cluster")
add_palette_to_uns(adata, "survivorship")
add_palette_to_uns(adata, "group", ["PT", "30HR", "250", "400"])
add_palette_to_uns(adata, "sample")

# %%
# TODO: should this be elsewhere or removed? does anything consume this file?
esam_expr = pl.from_pandas(
    sc.get.obs_df(
        adata[
            (adata.obs["sample"] == "Pretreatment") & (~adata.obs["barcode"].isna())
        ].copy(),
        keys=["barcode", "clone-cluster", "ESAM"],
    )
)

esam_expr.write_csv(project.paths.data / "pretreatment-esam-expression-by-barcode.csv")

# %%
hp.get_obs(adata).write_csv(project.paths.data / "final-obs.csv")
adata.write(dataset.final)

# %%
# ---
