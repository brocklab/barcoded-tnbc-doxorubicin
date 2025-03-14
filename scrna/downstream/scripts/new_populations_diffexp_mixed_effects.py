import argparse

import helpers as hp
import polars as pl
import scanpy as sc
from helpers.scrna.diffexp import find_de_MAST_barcode


def subset_anndata(adata_, post_sample, bgID):
    """
    return anndata containing pretreatment and
    one-postreatment sample
    """

    clone_cluster = (
        hp.get_obs(adata_)
        .filter(pl.col("bgID") == bgID)
        .select("clone-cluster")
        .unique()
        .item()
    )
    in_pretreatment_clone_cluster = (pl.col("sample") == "Pretreatment") & (
        pl.col("clone-cluster") == clone_cluster
    )
    in_post_sample_barcode = (pl.col("bgID") == bgID) & (
        pl.col("sample") == post_sample
    )

    cells = (
        hp.get_obs(adata_)
        .filter(in_pretreatment_clone_cluster | in_post_sample_barcode)["cell"]
        .to_list()
    )

    asub = adata_[cells, :].copy()
    return asub


def out_file(sample, bgID):
    return f"deg-new-population-{sample}-bgID-{bgID}-vs-Pretreatment.csv"


def main():
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

    project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
    dataset = hp.data.get_sc_dataset("231-1KB3")
    adata = sc.read(dataset.final)

    assert bgID in adata.obs["bgID"].cat.categories, f"bgID {bgID} not found in adata"
    assert sample in adata.obs["sample"].cat.categories, (
        f"sample {sample} not found in adata"
    )
    adata_ = subset_anndata(adata, sample, bgID)

    assert adata_.shape[0] > 0, "no cells passed filtering"
    (project.paths.data / "degs").mkdir(exist_ok=True)
    project.log.info(
        f"performing differential expression for {sample=}/{bgID=},"
        f"with anndata:\n {adata_}"
    )
    project.log.info("running differential expression")
    diffexp = find_de_MAST_barcode(adata_, "sample", "Pretreatment")
    diffexp.write_csv(project.paths.data / "degs" / out_file(sample, bgID))


# %%
if __name__ == "__main__":
    main()
# %%
# ---
# %%
