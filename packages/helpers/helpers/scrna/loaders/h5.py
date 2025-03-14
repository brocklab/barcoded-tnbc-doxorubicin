import argparse
import sys
from pathlib import Path

import scanpy as sc
from scipy.sparse import csr_matrix

import helpers as hp

hp.set_project("scrna/downstream")


class H5DataLoader:
    def __init__(self, h5: Path):
        self.h5 = h5

    def load(self):
        adata = sc.read_10x_h5(self.h5)
        adata.var_names_make_unique()
        # save space with sparse matrix
        adata.X = csr_matrix(adata.X)
        return adata


def main():
    parser = argparse.ArgumentParser()
    logger = hp.get_logger("h5-loader")

    for flag, help in (
        ("--h5", "path to h5 file"),
        ("--output", "path to write raw anndata"),
    ):
        parser.add_argument(flag, help=help, required=True, type=Path)

    parser.add_argument(
        "--sample-list", help="comma seperated list of samples (order matters)"
    )
    parser.add_argument("--samples", help="comma seperated list of samples to keep")

    parser.add_argument("dataset", help="name of dataset")

    args = parser.parse_args()

    if not args.h5.is_file():
        logger.critical(f"error: {args.h5} does not exist")
        sys.exit(1)

    sample_map = {i: sample for i, sample in enumerate(args.sample_list.split(","), 1)}

    logger.info(
        f"""
Arguments:
  H5: {args.h5}
  Output h5ad: {args.output}
  Sample Mapping: {sample_map}
  Samples: {args.samples}
"""
    )

    logger.info("reading raw data")
    adata = H5DataLoader(args.h5).load()

    logger.info("mapping samples to cell barcodes")
    if adata.obs.index.str.split("-").str[1].nunique() != len(sample_map):
        logger.error(
            "Expected anndata to and sample list to have the same number of samples"
        )
        sys.exit(1)

    adata.obs["sample"] = (
        adata.obs.index.str.split("-").str[1].astype(int).map(sample_map)
    )
    logger.info(f"subsetting anndata on samples: {args.samples}")

    adata = adata[adata.obs["sample"].isin(args.samples.split(","))].copy()
    if adata.shape[0] == 0:
        logger.error("anndata has no remaining observations")
        sys.exit(1)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"writing anndata to {args.output}")
    adata.write(args.output)


if __name__ == "__main__":
    main()
