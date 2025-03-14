# %%
import argparse

import helpers as hp
import polars as pl
import scanpy as sc
from helpers.scrna.diffexp import find_de_MAST_barcode


def out_file(sample, cluster, random):
    rand = "-randomized" if random else ""
    return f"deg-{sample}-cluster-{cluster}-survivorship-low-vs-high{rand}.csv"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", help="sample name", required=True, type=str)
    parser.add_argument("--cluster", help="cluster number", required=True, type=str)
    parser.add_argument("--random", help="with randomized labels", action="store_true")
    args = parser.parse_args()
    sample, cluster, random = args.sample, args.cluster, args.random

    project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
    dataset = hp.data.get_sc_dataset("231-1KB3")
    adata = sc.read(dataset.final)

    if random:
        project.log.info("including randomization!")
        randomized_survivorship = pl.read_csv(
            project.paths.data / "randomized-barcode-survivorship.csv"
        )
        adata.obs["survivorship"] = adata.obs["barcode"].map(
            dict(randomized_survivorship.iter_rows())
        )

    assert cluster in adata.obs["clone-cluster"].cat.categories, (
        f"cluster {cluster} not found in adata"
    )
    assert sample in adata.obs["sample"].cat.categories, (
        f"sample {sample} not found in adata"
    )

    default_filters = adata.obs["clone-specific"]
    adata_ = adata[
        default_filters
        & (adata.obs["clone-cluster"] == cluster)
        & (adata.obs["sample"] == sample)
    ].copy()
    assert adata_.shape[0] > 0, "no cells passed filtering"
    (project.paths.data / "degs").mkdir(exist_ok=True)
    project.log.info(
        f"performing differential expression for {cluster=}/{sample=},"
        f"with anndata:\n {adata_}"
    )
    project.log.info("running differential expression")
    diffexp = find_de_MAST_barcode(
        adata_,
        "survivorship",
        "low",
    )
    diffexp.write_csv(project.paths.data / "degs" / out_file(sample, cluster, random))


# %%
if __name__ == "__main__":
    main()
# %%
# ---
