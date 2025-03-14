import scanpy as sc

import helpers as hp

from . import add_var_categories, is_outlier, plot_qc_plots
from .doublets import run_scDblFinder

project = hp.Project("scrna/downstream", "scrna.qc", outs="outs/qc")
project.setup_scanpy()


def main():
    args = hp.scrna.dataset_parser().parse_args()
    dataset = hp.data.get_sc_dataset(args.dataset)

    project.log.info(f"Performing QC for dataset: {dataset}")

    adata = sc.read(dataset.raw)

    project.log.info(adata)

    add_var_categories(adata)

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    plot_qc_plots(adata, dataset.name)

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    project.log.info(f"Outlier cells by counts: \n{adata.obs.outlier.value_counts()}")
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 20
    )
    project.log.info(f"Outlier cells by mt: \n{adata.obs.mt_outlier.value_counts()}")

    project.log.info(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    project.log.info(
        f"Number of cells after filtering of low quality cells: {adata.n_obs}"
    )

    plot_qc_plots(adata, dataset.name, suffix="_post_qc")

    run_scDblFinder(adata)

    project.log.info(f"Number of genes prior to filtering: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=20)
    project.log.info(f"Number of genes after filtering: {adata.n_vars}")

    project.log.info(f"saving anndata to {dataset.postqc}")

    adata.write(dataset.postqc)


if __name__ == "__main__":
    main()
