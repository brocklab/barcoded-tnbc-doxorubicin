from pathlib import Path

import altair as alt
import matplotlib.pyplot as plt
import polars as pl
import scanpy as sc
from scipy.cluster.hierarchy import (
    dendrogram,
    fcluster,
    linkage,
    set_link_color_palette,
)

import helpers as hp

project = hp.Project("scrna/downstream", "scrna.barcodes", outs="outs/barcodes")
hp.setup_altair()
project.setup_scanpy()


@project.save_plot()
def clonal_cluster_per_sample(obs):
    return (
        alt.Chart(obs.group_by("clone-cluster", "sample").agg(pl.len().alias("cells")))
        .mark_bar()
        .encode(
            alt.X("sample"),
            alt.Y("cells"),
            color=(alt.Color("clone-cluster").title("clonal cluster")),
        )
        .properties(width=200)
    )


def plot_barcode_dendrogram(Z, labels, threshold):
    # TODO: improve colors
    fig = plt.figure()
    set_link_color_palette(["#9a46ab", "#127030"])
    dendrogram(
        Z,
        truncate_mode="level",
        p=4,
        labels=labels,
        orientation="right",
        color_threshold=threshold,
        above_threshold_color="#bcbddc",
    )
    fig.savefig(
        project.paths.outs / "dendrogram-barcode-cluster.svg", bbox_inches="tight"
    )

    fig.show()


def generate_barcode_linkage(adata):
    agg = sc.get.aggregate(
        adata[(~adata.obs["barcode"].isna()) & (adata.obs["sample"] == "Pretreatment")],
        by="barcode",
        func=["mean"],
    )

    X = agg.layers["mean"]
    labels = agg.obs["barcode"].to_list()
    Z = linkage(X, "ward")
    return Z, labels


def assign_clonal_cluster(adata, Z, labels, t):
    label_cluster_map = dict(
        zip(
            labels,
            map(str, map(int, fcluster(Z, t=t, criterion="distance"))),
        )
    )
    adata.obs["clone-cluster"] = (
        adata.obs["barcode"].map(label_cluster_map).fillna("unknown").astype("category")
    )


def surviving_barcodes():
    treated = hp.data.targeted_final.load()
    return treated.filter(
        (pl.col("sample") != "PT-0") & (pl.col("barcode") != "unknown")
    )["barcode"].to_list()


def main():
    project.log.info("BEGIN: single cell barcode clustering and assignment.")
    parser = hp.scrna.dataset_parser()
    parser.add_argument(
        "--assignments",
        help="cell barcode assignments",
        type=Path,
        default=hp.data.cell_barcodes_final.path,
    )
    parser.add_argument(
        "--threshold", help="barcode cluster threshold cutoff", type=int, default=200
    )
    args = parser.parse_args()

    dataset = hp.data.get_sc_dataset(args.dataset)
    adata = sc.read_h5ad(
        dataset.clustered,
    )

    adata.obs = (
        adata.obs.join(
            pl.read_csv(
                args.assignments,
            )
            .drop("sample")
            .to_pandas()
            .set_index("cell"),
            how="left",
        )
        .assign(
            survivor=lambda df: df["barcode"]
            .isin(surviving_barcodes())
            .map({True: "yes", False: "no"})
            .astype("category")
        )
        .assign(
            **{
                col: (
                    lambda x: lambda df: df[x]
                    .map({True: "yes", False: "no"})
                    .astype("category")
                )(col)
                for col in ("multiple-integration", "barcode-doublet", "cell-doublet")
            }
        )
    )

    project.log.info(
        "Removing {cells} cells that are perceived doublets from barcodes".format(
            cells=(adata.obs["barcode-doublet"] == "yes").sum()
        )
    )
    adata = adata[~(adata.obs["barcode-doublet"] == "yes")].copy()

    Z, labels = generate_barcode_linkage(adata)
    plot_barcode_dendrogram(Z, labels, args.threshold)

    assign_clonal_cluster(adata, Z, labels, args.threshold)

    project.log.info(adata.obs.groupby("sample").agg({"clone-cluster": "value_counts"}))

    clonal_cluster_per_sample(hp.get_obs(adata))

    project.log.info(f"writing h5ad to {dataset.barcoded}")
    adata.write(dataset.barcoded)
    project.log.info("END: single cell barcode clustering and assignment.")


if __name__ == "__main__":
    main()
