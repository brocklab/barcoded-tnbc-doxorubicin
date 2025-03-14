import scanpy as sc

import helpers as hp

from . import normalize_counts, plot_normalized_counts

project = hp.Project("scrna/downstream", "scrna.normalize", outs="outs/normalize")
project.setup_scanpy()


def main():
    args = hp.scrna.dataset_parser().parse_args()
    dataset = hp.data.get_sc_dataset(args.dataset)

    adata = sc.read(dataset.postqc)

    project.log.info("Removing doublets")
    adata = adata[adata.obs["scDblFinder_class"] == "singlet"].copy()

    project.log.info("saving raw counts in counts layer")
    adata.layers["counts"] = adata.X.copy()

    project.log.info(f"Performing normalization for dataset: {dataset.name}")
    normalize_counts(adata)

    plot_normalized_counts(adata, dataset.name)

    project.log.info("updating counts with normalized counts")
    # adata.X = adata.layers.pop("log1p_norm")
    adata.X = adata.layers["log1p_norm"]

    project.log.info("setting raw annadata")
    adata.raw = adata.copy()

    project.log.info(f"Writing anndata to {dataset.normalized}")

    adata.write(dataset.normalized)


if __name__ == "__main__":
    main()
