import altair as alt
import helpers as hp
import polars as pl
import scanpy as sc

hp.setup_altair()
project = hp.Project("scrna/barcodes", "scrna.barcodes.filter")

dataset = hp.data


@project.save_plot()
def filtered_cells_bar(df):
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            alt.X("sample").axis(ticks=False),
            alt.Y("# of cells"),
            color=alt.Color("keep").scale(range=("red", "blue")),
        )
    ).properties(width=150)


def main():
    dataset = hp.data.get_sc_dataset("231-1KB3")
    project.log.info("BEGIN: filtering barcodes.")
    samples = dataset.samples
    sample_idx = pl.DataFrame(
        ((i, sample.id) for i, sample in enumerate(samples, 1)),
        schema=("sample-idx", "sample"),
    )
    barcodes = (
        pl.read_csv(project.paths.data / "collapsed-barcodes.csv")
        .drop("barcode")
        .rename({"collapsed-barcode": "barcode"})
    )

    adata = sc.read(dataset.postqc)

    obs = pl.from_pandas(adata.obs.reset_index()).select(
        pl.col("sample").cast(str), pl.col("index").alias("cell")
    )
    barcodes = (
        barcodes.join(sample_idx, on="sample")
        .with_columns(
            (cell := pl.col("cell") + "-" + pl.col("sample-idx").cast(str)),
            cell.is_in(obs["cell"].to_list()).alias("keep"),
        )
        .drop("sample-idx")
    )
    project.log.info(barcodes.group_by("sample", "keep").agg(pl.col("cell").n_unique()))

    filtered_cells_bar(
        barcodes.group_by("sample", "keep").agg(
            pl.col("cell").n_unique().alias("# of cells")
        )
    )
    barcodes = barcodes.filter(pl.col("keep")).drop("keep")

    project.log.info(
        f"writing two csvs to {project.paths.data}/\n"
        "- collapsed-filtered-barcodes.csv\n- barcodes-pre-clustering.tsv"
    )

    barcodes.write_csv(project.paths.data / "collapsed-filtered-barcodes.csv")
    barcodes.select("barcode").write_csv(
        project.paths.data / "barcodes-pre-clustering.tsv", include_header=False
    )
    project.log.info("END: filtering barcodes.")


if __name__ == "__main__":
    main()
