import scanpy as sc

import helpers as hp

RANDOM_STATE = 231


def get_highly_variable_genes(adata, dataset):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, save=f"_{dataset}.svg")


def compute_dimensionality_reductions(
    adata,
):
    sc.pp.pca(
        adata,
        svd_solver="arpack",
        mask_var="highly_variable",
        random_state=RANDOM_STATE,
    )
    sc.pp.neighbors(adata, random_state=RANDOM_STATE)
    sc.tl.umap(adata, random_state=RANDOM_STATE)


def plot_dimensionality_reductions(adata, dataset):
    sc.pl.pca_scatter(adata, color="total_counts", save=f"_counts_{dataset}.svg")
    sc.pl.umap(
        adata,
        color=["total_counts", "pct_counts_mt"],
        save=f"_counts_pct_mt_{dataset}.png",
    )


def main():
    project = hp.Project("scrna/downstream", "scrna.cluster", outs="outs/cluster")
    project.setup_scanpy()

    parser = hp.scrna.dataset_parser()
    parser.add_argument(
        "--resolution",
        help="leiden resolution",
        type=float,
        default=0.1,  # low resolution default designed for 231-1KB3
    )
    args = parser.parse_args()
    dataset = hp.data.get_sc_dataset(args.dataset)

    adata = sc.read(dataset.normalized)
    project.log.info("prior to clustering, computing phase score")
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=hp.data.cell_cycle_genes.s_genes,
        g2m_genes=hp.data.cell_cycle_genes.g2m_genes,
    )

    project.log.info("computing highly variable genes")
    get_highly_variable_genes(adata, dataset.name)

    project.log.info("computing dimensionality reductions")
    compute_dimensionality_reductions(adata)

    project.log.info("plotting dimensionality reductions")
    plot_dimensionality_reductions(adata, dataset.name)

    project.log.info("computing leiden clusters")
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        # future scanpy defaults
        flavor="igraph",
        n_iterations=2,
    )

    project.log.info(f"Writing anndata to {dataset.clustered}")
    adata.write(dataset.clustered)


if __name__ == "__main__":
    main()
