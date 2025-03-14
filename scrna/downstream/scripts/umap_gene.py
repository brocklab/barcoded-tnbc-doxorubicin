# %%
import argparse

import helpers as hp
import scanpy as sc
from helpers.plotting.umap import Umap


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", help="gene name", required=True, type=str)
    args = parser.parse_args()
    gene = args.gene

    project = hp.Project("scrna/downstream", "scrna.downstream.umap", outs="outs/umaps")
    dataset = hp.data.get_sc_dataset("231-1KB3")
    project.log.info("loading data")
    adata = sc.read(dataset.final)

    project.log.info(f"plotting umap for: {gene}")
    Umap(adata, key=gene).plot(save=project.paths.outs / f"umap-{gene}.png")


# %%
if __name__ == "__main__":
    main()
# %%
# ---
