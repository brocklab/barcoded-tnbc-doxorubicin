import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from helpers.utils import get_logger

TYPES = {"str": str, "int32": np.int32, "float64": np.float64}


class RawDataLoader:
    def __init__(
        self, umis: Path, metadata: Path, metadata_kws: Dict[str, Any] | None = None
    ):
        self.umis = umis
        self.metadata = metadata
        self.metadata_kws = metadata_kws

    def _load_metadata(self):
        return pd.read_csv(
            self.metadata, **(self.metadata_kws if self.metadata_kws else {})
        )

    def load(self, transpose: bool = False):
        adata = anndata.read_text(self.umis)
        metadata = self._load_metadata()

        # save space with sparse matrix
        adata.X = csr_matrix(adata.X)

        # determine if transpose is necessary using metadata size?
        if transpose:
            adata = adata.T

        adata.obs = adata.obs.join(metadata)
        return adata


def main():
    parser = argparse.ArgumentParser()
    logger = get_logger("raw-loader")

    for flag, help in (
        ("--umis", "path to umi counts as txt file"),
        ("--metadata", "path to metdata as txt file"),
        ("--output", "path to write raw anndata"),
    ):
        parser.add_argument(flag, help=help, required=True, type=Path)
    parser.add_argument(
        "--metadata_kws", help="json string with kw arguments", type=str, default="{}"
    )
    parser.add_argument(
        "--filter-na", help="column that must have a value in adata.obs", type=str
    )
    parser.add_argument("--transpose", help="transpose umis", action="store_true")
    args = parser.parse_args()

    for p in (args.umis, args.metadata):
        if not p.is_file():
            logger.critical(f"error: {p} does not exist")
            sys.exit(1)

    kwargs = json.loads(args.metadata_kws)

    if "dtype" in kwargs:
        kwargs["dtype"] = {k: TYPES[v] for k, v in kwargs["dtype"].items()}

    logger.info(
        f"""
Arguments:
  Umis: {args.umis}
  Transpose: {args.transpose}
  Metadata: {args.metadata}
  Metdata Loading Args: {kwargs}
  Output h5ad: {args.output}
"""
    )

    logger.info("reading raw data")

    adata = RawDataLoader(args.umis, args.metadata, kwargs).load(args.transpose)
    args.output.parent.mkdir(exist_ok=True, parents=True)

    if args.filter_na:
        logger.info(f"filtering on {args.filter_na}")
        logger.info(adata)
        adata = adata[~adata.obs[args.filter_na].isna()].copy()
        logger.info(adata)

    logger.info(f"writing anndata to {args.output}")

    adata.write(args.output)


if __name__ == "__main__":
    main()
