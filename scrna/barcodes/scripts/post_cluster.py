import helpers as hp
import polars as pl

project = hp.Project("scrna/barcodes", "scrna.barcodes.post_cluster")


def main():
    project.log.info("BEGIN: post-cluster processing barcodes.")
    post_cluster_barcodes = (
        pl.read_csv(
            project.paths.data / "barcodes-post-clustering.tsv",
            separator="\t",
            has_header=False,
            new_columns=("centroid", "cluster"),
            columns=[0, 2],
        )
        .with_columns(pl.col("cluster").str.split(","))
        .explode("cluster")
        .rename({"cluster": "barcode"})
    )

    barcodes = pl.read_csv(project.paths.data / "collapsed-filtered-barcodes.csv").join(
        post_cluster_barcodes, on="barcode", how="full"
    )
    project.log.info("total centroids vs unique barcodes:")
    project.log.info(
        barcodes.group_by("sample").agg(
            pl.col("barcode").n_unique(), pl.col("centroid").n_unique()
        )
    )
    project.log.info(
        f"writing csv to {project.paths.data / 'collapsed-filtered-clustered-barcodes.csv'}"
    )
    (
        barcodes.drop("barcode")
        .rename({"centroid": "barcode"})
        .write_csv(project.paths.data / "collapsed-filtered-clustered-barcodes.csv")
    )
    project.log.info("END: post-cluster processing barcodes.")


if __name__ == "__main__":
    main()
