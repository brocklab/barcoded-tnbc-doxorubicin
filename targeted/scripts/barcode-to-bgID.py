# %%
import helpers as hp
import polars as pl


def assign_barcode_ids(df):
    "assign a unique ID to each barcode bgXXX"
    rank_col = (
        pl.col("percent")
        .rank(method="ordinal", descending=True)
        .alias("rank")
        .cast(int)
    )
    ids = (
        df.sort("percent", "barcode")
        .with_columns(("bg" + (rank_col).cast(str).str.zfill(3)).alias("bgID"))
        .select("barcode", "bgID")
    )

    assert ids["barcode"].n_unique() == ids["bgID"].n_unique(), (
        "barcode id assignment has unexpected overlaps"
    )

    return df.with_columns(
        pl.col("barcode")
        .replace_strict(old=ids["barcode"], new=ids["bgID"])
        .alias("bgID")
    )


project = hp.Project("targeted")
library = pl.read_csv(project.paths.data / "counts-library.tsv", separator="\t")

assign_barcode_ids(library).select("barcode", "bgID").write_csv(
    project.paths.data / "barcode-bgID.csv"
)
