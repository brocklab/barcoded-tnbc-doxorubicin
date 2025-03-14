# %%
import helpers as hp
import polars as pl
import scanpy as sc

# %%
project = hp.Project("scrna/downstream")
hp.setup_altair()
project.setup_scanpy()
# %%
adata = sc.read(hp.data.get_sc_dataset("231-1KB3").final)
obs = hp.get_obs(adata)


# %%


def with_totals(df: pl.DataFrame, col="sample") -> pl.DataFrame:
    """
    given some single observation column return a dataframe
    containing totals of all other columns with entry in category as "total"
    Args:
        df: a dataframe with a single categorical column
        col: name of categorical column
    Returns;
        dataframe with all but "col" sum()'ed

    """
    enum_dtype = pl.Enum(df["sample"].unique().cast(str).append(pl.Series(["totals"])))
    return pl.concat(
        [
            df.with_columns(pl.col(col).cast(enum_dtype)),
            df.select(pl.all().exclude(col))
            .sum()
            .with_columns(pl.lit("totals").cast(enum_dtype).alias(col)),
        ],
        how="diagonal",
    ).sort(col)


# TODO: add to helpers?
class Report:
    data: str = ""

    # def __init__(self) -> None:
    def add(self, *data: str) -> None:
        self.data += " ".join(map(repr, data))

    def add_line(self, *data: str) -> None:
        print(*data)
        self.add(*data)
        self.add("\n")

    def add_header(self, header: str, level=1) -> None:
        print("#" * level, header)
        self.add("#" * level + " " + header + "\n\n")

    def add_table(self, data: pl.DataFrame) -> None:
        print(data)
        if "sample" in data.columns:
            data = data.sort("sample")
        with pl.Config(
            tbl_formatting="MARKDOWN",
            tbl_hide_dataframe_shape=True,
            tbl_hide_column_data_types=True,
        ):
            self.add("\n" + repr(data) + "\n")

    def show(self) -> None:
        print(self.data)


report = Report()
# %%
report.add_header("Cells + Reads")
# %%
report.add_line("per sample:")
report.add_table(
    obs.group_by("sample").agg(
        pl.len().alias("cells (N)"),
        pl.col("total_counts").median().cast(int),
        pl.col("n_genes_by_counts").median().cast(int),
    )
)
# %%
report.add_line("totals:")
report.add_table(
    obs.select(
        pl.len().alias("cells (N)"),
        pl.col("total_counts").median().cast(int),
        pl.col("n_genes_by_counts").median().cast(int),
    )
)
# %%
# %%
report.add_header("Barcodes")
# %%
report.add_line("percent of cells per sample by number of assigned barcodes:")
report.add_table(
    (
        obs.group_by("sample", "n_barcodes")
        .len()
        .select(
            "sample",
            pl.col("n_barcodes").fill_null(0).cast(int),
            ((pl.col("len") / pl.col("len").sum().over("sample")) * 100)
            .alias("percent")
            .round(3),
        )
    )
    .pivot(index="sample", on="n_barcodes", sort_columns=True)
    .sort("sample")
    .fill_null(0)
)
# %%
report.show()
# %%
report.add_line("unique barcode (including multiple barcodes) per sample:")
report.add_table(
    obs.filter(pl.col("barcode").is_not_null())
    .group_by("sample")
    .agg(pl.col("barcode").n_unique().alias("unique barcodes"))
    .sort("sample")
)

# %%
report.add_line("unique barcodes (in only single barcode cells) per sample:")
report.add_table(
    obs.filter(pl.col("n_barcodes") == 1)
    .group_by("sample")
    .agg(pl.col("barcode").n_unique().alias("unique barcodes"))
    .sort("sample")
)

# %%
report.add_line(
    "total percent cells assigned single unique barcode:",
    obs.select(((pl.col("n_barcodes") == 1).sum() / pl.len() * 100).round(2)).item(),
)
# %%
report.show()
# %%
# %%
