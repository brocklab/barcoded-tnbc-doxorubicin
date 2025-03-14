# %%
import argparse

import helpers as hp
import polars as pl
from helpers.scrna.gsea import run_gsea_prerank


def main():
    project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster", help="cluster number", required=True, type=str)
    args = parser.parse_args()
    cluster = args.cluster

    degs = pl.read_csv(
        project.paths.data
        / "degs"
        / f"deg-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
    )

    gsea = run_gsea_prerank(degs)
    (project.paths.data / "gsea").mkdir(exist_ok=True)
    gsea.write_csv(
        project.paths.data
        / "gsea"
        / f"gsea-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
    )


# %%
if __name__ == "__main__":
    main()

# %%
# ---
