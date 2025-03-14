# %%
import argparse

import helpers as hp
import scanpy as sc
from helpers.scrna.diffexp import find_de_MAST_barcode


def out_file(sample):
    return f"deg-{sample}-clone-cluster-1-vs-2.csv"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", help="sample name", required=True, type=str)
    args = parser.parse_args()
    sample = args.sample

    project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
    dataset = hp.data.get_sc_dataset("231-1KB3")
    adata = sc.read(dataset.final)

    assert sample in adata.obs["sample"].cat.categories, (
        f"sample {sample} not found in adata"
    )

    default_filters = adata.obs["clone-specific"] & (
        adata.obs["clone-cluster"] != "unknown"
    )
    adata_ = adata[default_filters & (adata.obs["sample"] == sample)].copy()

    assert adata_.shape[0] > 0, "no cells passed filtering"
    (project.paths.data / "degs").mkdir(exist_ok=True)
    project.log.info(
        f"performing differential expression for {sample=},with anndata:\n {adata_}"
    )
    project.log.info("running differential expression")
    diffexp = find_de_MAST_barcode(
        adata_,
        "clone-cluster",
        "1",
    )
    diffexp.write_csv(project.paths.data / "degs" / out_file(sample))


# %%
if __name__ == "__main__":
    main()
# %%
# ---
