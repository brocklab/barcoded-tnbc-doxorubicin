# %%
import argparse

import helpers as hp
import scanpy as sc
from helpers.scrna.diffexp import find_de_MAST_barcode


def out_file(cluster):
    return f"deg-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster", help="cluster number", required=True, type=str)
    args = parser.parse_args()
    cluster = args.cluster

    project = hp.Project("scrna/downstream", "scrna.downstream.barcodes")
    dataset = hp.data.get_sc_dataset("231-1KB3")
    adata = sc.read(dataset.final)
    assert cluster in adata.obs["clone-cluster"].cat.categories, (
        f"cluster {cluster} not found in adata"
    )

    default_filters = adata.obs["clone-specific"] & adata.obs["sample"].isin(
        ("Pretreatment", "Dox-550-30HR")
    )
    adata_ = adata[default_filters & (adata.obs["clone-cluster"] == cluster)].copy()

    assert adata_.shape[0] > 0, "no cells passed filtering"
    (project.paths.data / "degs").mkdir(exist_ok=True)
    project.log.info(
        f"performing differential expression for {cluster=},with anndata:\n {adata_}"
    )
    project.log.info("running differential expression")
    diffexp = find_de_MAST_barcode(
        adata_,
        "sample",
        "Pretreatment",
    )
    diffexp.write_csv(project.paths.data / "degs" / out_file(cluster))


# %%
if __name__ == "__main__":
    main()
# %%
# ---
