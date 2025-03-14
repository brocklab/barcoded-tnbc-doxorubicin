# %%
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted", "targeted.treated")
hp.setup_altair()
# %%
df = hp.data.targeted_final.load()
post_barcode_df = df.filter(pl.col("barcode") != "unknown", pl.col("sample") != "PT-0")
# %%
df.head()


# %%
# add to helpers?
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


# %%
report = Report()
# %%
report.add_header("Barcode Metrics")
report.add_line(
    "total unique barcodes at PT:",
    df.filter(pl.col("sample") == "PT-0")["barcode"].n_unique(),
)
# %%
# number of samples where top 3 clones >= 90% of population

n_three_clone_samples = (
    post_barcode_df.sort("percent", descending=True)
    .group_by("sample", maintain_order=True)
    .head(3)
    .group_by("sample")
    .agg(pl.col("percent").sum() > 90)["percent"]
    .sum()
    # .group_by("sample")
)
report.add_line(
    f"{n_three_clone_samples}/{df.filter(pl.col('sample') != 'PT-0')['sample'].n_unique()}"
    "samples are >= 90% comprised of just 3 clones"
)

# %%
report.add_line(
    "number of clones in more than one sample:",
    (
        post_barcode_df.group_by("barcode")
        .agg(pl.col("sample").n_unique() > 1)["sample"]
        .sum()
    ),
)
# %%
report.add_line(
    "number of clones at >90%, in at least one sample:",
    post_barcode_df.filter(pl.col("percent") > 90)["barcode"].n_unique(),
)
# %%
report.show()
# %%
