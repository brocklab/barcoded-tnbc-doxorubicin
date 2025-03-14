# %%
import argparse

import helpers as hp
import polars as pl
from helpers.scrna.gsea import run_gsea_prerank


def post_process_gsea(gsea):
    split_term = pl.col("full term").str.split(" (")
    return gsea.rename({"Term": "full term"}).with_columns(
        split_term.list.get(-1).alias("GO ID"), split_term.list.get(0).alias("term")
    )


def main():
    project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", help="sample name", required=True, type=str)
    parser.add_argument("--bgID", help="barcode id", required=True, type=str)
    args = parser.parse_args()
    (
        sample,
        bgID,
    ) = (
        args.sample,
        args.bgID,
    )

    degs = pl.read_csv(
        project.paths.data
        / "degs"
        / f"deg-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    )

    gsea = run_gsea_prerank(degs)
    (project.paths.data / "gsea").mkdir(exist_ok=True)
    gsea.write_csv(
        project.paths.data
        / "gsea"
        / f"gsea-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    )


# %%
if __name__ == "__main__":
    main()

# %%
# ---
