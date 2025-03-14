from collections import namedtuple
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import List

import anndata
import anndata as ad
import polars as pl
import scanpy as sc
from pandas.util import hash_pandas_object
from scipy.sparse import csr_matrix

# HARDCODED
ROOT_PATH = Path(__file__).parent.parent.parent.parent.parent


@dataclass
class Sample:
    id: str
    batch: str
    name: str
    group: str


@dataclass
class SingleCellDataSet:
    name: str
    samples: List[Sample]
    counts_10X_path: Path

    def __post_init__(self):
        self._data_dir = ROOT_PATH / "scrna" / "downstream" / "data"

    # @property
    # def h5(self) -> Path:
    #     return self._data_dir / "raw" / self.name / "filtered_feature_bc_matrix.h5"

    def _filename(self, property) -> Path:
        return self._data_dir / "h5ads" / f"{self.name}_{property}.h5ad"

    def _sample_10X_mtx_path(self, sample: Sample) -> Path:
        if sample not in self.samples:
            raise ValueError(
                f"Unknown sample: {sample}, expected one of: "
                + ",".join(map(str, self.samples))
            )
        return self.counts_10X_path / sample.id / "outs" / "filtered_feature_bc_matrix"

    def _load_sample(self, sample: Sample, idx: int) -> ad.AnnData:
        adata = sc.read_10x_mtx(self._sample_10X_mtx_path(sample))
        adata.var_names_make_unique()
        adata.obs.index = adata.obs.index.str.split("-").str[0] + "-" + str(idx)
        adata.obs["sample-id"] = sample.id
        adata.obs["sample"] = sample.name
        adata.obs["batch"] = sample.batch
        adata.obs["group"] = sample.group

        # save space with sparse matrix
        adata.X = csr_matrix(adata.X)
        return adata

    def load_raw(self) -> ad.AnnData:
        adatas = [
            self._load_sample(sample, i) for i, sample in enumerate(self.samples, 1)
        ]
        return ad.concat(
            adatas,
            join="outer",
        )

    @property
    def sample_ids(self) -> List[str]:
        return [sample.id for sample in self.samples]

    @property
    def sample_10X_mtx_paths(self) -> List[Path]:
        return [
            self._sample_10X_mtx_path(sample) / f
            for sample in self.samples
            for f in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
        ]

    @property
    def raw(self) -> Path:
        return self._filename("raw")

    @property
    def postqc(self) -> Path:
        return self._filename("postqc")

    @property
    def normalized(self) -> Path:
        return self._filename("postqc-normalized")

    @property
    def clustered(self) -> Path:
        return self._filename("postqc-normalized-clustered")

    @property
    def barcoded(self) -> Path:
        return self._filename("postqc-normalized-clustered-barcoded")

    @property
    def final(self) -> Path:
        return self._filename("postqc-normalized-clustered-barcoded-final")


_sc_datasets = {
    "231-1KB3": SingleCellDataSet(
        name="231-1KB3",
        samples=[
            # batch 1
            Sample("PT-1", "1", "Pretreatment", "PT"),
            Sample("Dox", "1", "Dox-550-30HR", "30HR"),
            Sample("A4", "1", "Dox-250-4", "250"),
            Sample("A6", "1", "Dox-250-6", "250"),
            # batch 2
            Sample("PT-2", "2", "Pretreatment", "PT"),
            Sample("A1", "2", "Dox-250-1", "250"),
            Sample("B2", "2", "Dox-400-2", "400"),
            Sample("B5", "2", "Dox-400-4", "400"),
            Sample("B6", "2", "Dox-400-5", "400"),
        ],
        counts_10X_path=ROOT_PATH / "scrna" / "fastq2h5ad" / "counts",
    )
}


def get_sc_dataset(name: str) -> SingleCellDataSet:
    if name not in _sc_datasets:
        raise ValueError(
            f"unknown dataset: {name}, must be one of: {','.join(_sc_datasets.keys())}"
        )
    return _sc_datasets[name]


def raw_data_paths() -> List[Path]:
    # NOTE: HARDCODED
    sc_fastqs = sorted(
        (ROOT_PATH / "scrna" / "fastq2h5ad" / "data" / "raw").glob("**/*.fastq.gz")
    )
    targeted_fastqs = sorted((ROOT_PATH / "targeted" / "data" / "raw").glob("*.fastq*"))
    misc = [
        (ROOT_PATH / "targeted" / "data" / "DM2103-bfp.csv"),
    ]

    return sc_fastqs + targeted_fastqs + misc


# HARDCODED!!
@dataclass
class CsvDataSet:
    path: Path

    def load(
        self,
    ):
        return pl.read_csv(self.path)


targeted_final = CsvDataSet(ROOT_PATH / "targeted" / "data" / "targeted-final.csv")

obs_final = CsvDataSet(
    # TODO: change path to obs-final.csv?
    ROOT_PATH / "scrna/downstream/data/final-obs.csv"
)

cell_barcodes_final = CsvDataSet(
    ROOT_PATH / "scrna/barcodes/data/cell-barcode-assignments.csv",
)

barcode_bgID = CsvDataSet(
    ROOT_PATH / "targeted/data/barcode-bgID.csv",
)


def load_pycashier_outs(
    path: Path, columns: List[str] = ["barcode", "count", "percent", "sample"]
) -> pl.DataFrame:
    split_col = pl.col("sample").str.split("-").list
    df = (
        pl.read_csv(path, separator="\t")
        .with_columns(
            split_col.get(0).alias("job"),
            split_col.slice(1).list.join("-").alias("sample"),
        )
        .sort("job", "sample")
    ).select(*columns)
    return df


# cell_cycle_genes from:
# https://github.com/satijalab/seurat/blob/master/data/cc.genes.updated.2019.rda

cell_cycle_genes = namedtuple("CellCycleGenes", "s_genes g2m_genes")(
    [
        "MCM5",
        "PCNA",
        "TYMS",
        "FEN1",
        "MCM7",
        "MCM4",
        "RRM1",
        "UNG",
        "GINS2",
        "MCM6",
        "CDCA7",
        "DTL",
        "PRIM1",
        "UHRF1",
        "CENPU",
        "HELLS",
        "RFC2",
        "POLR1B",
        "NASP",
        "RAD51AP1",
        "GMNN",
        "WDR76",
        "SLBP",
        "CCNE2",
        "UBR7",
        "POLD3",
        "MSH2",
        "ATAD2",
        "RAD51",
        "RRM2",
        "CDC45",
        "CDC6",
        "EXO1",
        "TIPIN",
        "DSCC1",
        "BLM",
        "CASP8AP2",
        "USP1",
        "CLSPN",
        "POLA1",
        "CHAF1B",
        "MRPL36",
        "E2F8",
    ],
    [
        "HMGB2",
        "CDK1",
        "NUSAP1",
        "UBE2C",
        "BIRC5",
        "TPX2",
        "TOP2A",
        "NDC80",
        "CKS2",
        "NUF2",
        "CKS1B",
        "MKI67",
        "TMPO",
        "CENPF",
        "TACC3",
        "PIMREG",
        "SMC4",
        "CCNB2",
        "CKAP2L",
        "CKAP2",
        "AURKB",
        "BUB1",
        "KIF11",
        "ANP32E",
        "TUBB4B",
        "GTSE1",
        "KIF20B",
        "HJURP",
        "CDCA3",
        "JPT1",
        "CDC20",
        "TTK",
        "CDC25C",
        "KIF2C",
        "RANGAP1",
        "NCAPD2",
        "DLGAP5",
        "CDCA2",
        "CDCA8",
        "ECT2",
        "KIF23",
        "HMMR",
        "AURKA",
        "PSRC1",
        "ANLN",
        "LBR",
        "CKAP5",
        "CENPE",
        "CTCF",
        "NEK2",
        "G2E3",
        "GAS2L3",
        "CBX5",
        "CENPA",
    ],
)


class WrappedObs:
    """basic wrapper to make adata.obs cacheable"""

    def __init__(self, obs):
        self.obs = obs

    def __eq__(self, other):
        return self.obs.equals(other)

    def __hash__(self):
        series_hash = hash_pandas_object(self.obs)
        return hash((tuple(series_hash.index), tuple(series_hash.values)))


@lru_cache(maxsize=5)
def cached_obs2pl(obs: WrappedObs) -> pl.DataFrame:
    return pl.from_pandas(obs.obs.reset_index().rename(columns={"index": "cell"}))


def get_obs(adata: anndata.AnnData) -> pl.DataFrame:
    return cached_obs2pl(WrappedObs(adata.obs))
