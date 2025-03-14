# %%
import argparse

import helpers as hp
import scanpy as sc
from helpers.scrna.diffexp import find_de_MAST


def subset_anndata(adata_, sample, bgID):
    filters = adata_.obs["sample"].isin(["Pretreatment", sample]) & (
        adata_.obs["bgID"] == bgID
    )
    return adata_[filters].copy()


def out_file(sample, bgID):
    return f"deg-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"


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
    diffexp = find_de_MAST(adata_, "sample", "Pretreatment")
    diffexp.write_csv(project.paths.data / "degs" / out_file(sample, bgID))


# %%
if __name__ == "__main__":
    main()
# %%
# ---
