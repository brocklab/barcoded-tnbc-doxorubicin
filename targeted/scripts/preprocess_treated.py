# %%
import altair as alt
import helpers as hp
import polars as pl

# %%
project = hp.Project("targeted")
hp.setup_altair()

project.log.info("BEGIN: preprocessing barcode treated data.")


# %%
def calc_percent_cpm_log2cpm(df):
    return df.with_columns(
        (
            (prop := pl.col("count") / pl.col("count").sum().over(pl.col("sample")))
            * 100
        ).alias("percent"),
        (cpm := (prop * 1e6).alias("cpm")),
        cpm.log(2).alias("log2cpm"),
    )


# %%
lib_bgIDs = pl.read_csv(project.paths.data / "barcode-bgID.csv")
df = hp.data.load_pycashier_outs(project.paths.data / "counts-treated.tsv").select(
    "barcode", "count", "sample"
)
project.log.info(df.head())


# %%
allow_list = set(lib_bgIDs["barcode"]).intersection(
    df.filter(pl.col("sample") == "1KB3")["barcode"]
)
in_pt = (
    df.with_columns(pl.col("barcode").is_in(allow_list).alias("in-pt"))
    .group_by("sample")
    .agg(pl.col("in-pt").value_counts())
    .explode("in-pt")
    .unnest("in-pt")
    .with_columns(
        (pl.col("count") / pl.col("count").sum().over(pl.col("sample"))).alias(
            "proportion"
        ),
        (~pl.col("in-pt").alias("filtered")),
    )
)


# %%
@project.save_plot()
def proportion_barcodes_in_pretreatment(in_pt):
    return (
        alt.Chart(in_pt)
        .mark_bar(stroke="black")
        .encode(
            alt.X("proportion").title("proportion of unique barcodes").axis(offset=5),
            alt.Y("sample").axis(ticks=False, offset=5, domain=False),
            color=alt.Color("filtered").scale(range=["white", "cyan"]).title("removed"),
        )
        .properties(width=300)
    )


proportion_barcodes_in_pretreatment(in_pt)
# %%
project.log.info("barcodes in allow list")
project.log.info(
    df.select(pl.col("barcode").unique().is_in(allow_list))["barcode"].value_counts()
)
df = df.filter(pl.col("barcode").is_in(allow_list))
project.log.info(f"writing csv to {project.paths.data / 'treated-filtered.csv'}")
df.write_csv(project.paths.data / "treated-filtered.csv")
# %%
# treat 99 BFP+ by flow as 100% BFP cells
bfp = pl.read_csv(project.paths.data / "DM2103-bfp.tsv").with_columns(
    pl.when(pl.col("bfp") == 99).then(pl.lit(100)).otherwise(pl.col("bfp")).alias("bfp")
)

# C4 has ~ 15K remaining reads so we will filter
# it out here with a minimum of 50K valid reads
df = df.filter(
    pl.col("sample").is_in(
        df.group_by("sample")
        .agg(pl.col("count").sum())
        .filter(pl.col("count") > 50_000)["sample"]
    )
)

# %%
unknowns = (
    df.group_by("sample")
    .agg(pl.col("count").sum())
    .filter(pl.col("count") > 50_000)
    .join(bfp.filter(~pl.col("bfp").is_in([0, 100])), on="sample")
    .with_columns((pl.col("count") / (pl.col("bfp") / 100)).alias("new_count"))
    .select(
        barcode=pl.lit("unknown"),
        count=(pl.col("new_count") - pl.col("count")).cast(pl.Int64).alias("count"),
        sample=pl.col("sample"),
    )
)

# %%
df_unknown = pl.concat([df, unknowns]).with_columns(
    (pl.col("count") / pl.col("count").sum().over("sample") * 100).alias("percent")
)
# %%
df.group_by("sample").agg(pl.col("barcode").n_unique()).plot.bar(
    x="barcode", y="sample"
)

# %%
# apply one final abundance filter
# after scaling for noise/unknown
# barcodes then recalculate total percent
# df_unknown = df_unknown.filter(pl.col("percent") > 0.01).with_columns(
#     (pl.col("count") / pl.col("count").sum().over("sample") * 100).alias("percent")
# )
# %%


def assign_sample_names(df):
    # some samples didn't recover so we will reasign their replicate id's as necessary
    # to push the non-surviving samples to the end
    # e.g. B3 -> B6
    to_replace = {
        "B4": "B3",
        "B5": "B4",
        "B6": "B5",
        "B3": "B6",
        "C2": "C1",
        "C3": "C2",
        "C5": "C3",
        "C6": "C4",
        "C4": "C5",
        "C1": "C6",
    }

    group_col = (
        pl.col("sample")
        .str.slice(0, length=1)
        .replace({"A": "250", "B": "400", "C": "550", "1": "PT"})
    ).alias("group")
    replicate_col = (
        pl.when(pl.col("sample") != "1KB3")
        .then(pl.col("sample").replace(to_replace).str.slice(-1))
        .otherwise(pl.lit("0"))
    ).alias("replicate")

    # add sample metadata
    df = df.with_columns(
        group_col,
        replicate_col,
        pl.col("sample").alias("exp-id"),
        (group_col + pl.lit("-") + replicate_col).alias("sample"),
    )

    return df


df_unknown_samples = assign_sample_names(df_unknown)
# %%


def assign_missing_ids(df, lib_bgIDs):
    maxID = lib_bgIDs["bgID"].str.replace("bg", "").cast(int).max()

    return (
        df.filter(
            ~pl.col("barcode").is_in(lib_bgIDs["barcode"]),
            pl.col("barcode") != "unknown",
        )
        .group_by("barcode")
        .agg(pl.col("percent").mean())
        .with_columns(
            (
                pl.col("percent").rank(method="ordinal", descending=True).cast(int)
                + maxID
            ).alias("rank")
        )
        .select("barcode", ("bg" + pl.col("rank").cast(str)).alias("bgID"))
    )


def assign_barcode_ids(df, lib_bgIDs):
    bgIDs = pl.concat(
        (
            lib_bgIDs,
            assign_missing_ids(df_unknown_samples, lib_bgIDs),
            pl.DataFrame(dict(barcode=["unknown"], bgID=["unknown"])),
        )
    )

    return df.with_columns(
        pl.col("barcode")
        .replace_strict(old=bgIDs["barcode"], new=bgIDs["bgID"])
        .alias("bgID")
    )


df_unknown_samples_ids = assign_barcode_ids(df_unknown_samples, lib_bgIDs)
# %%
project.log.info(f"writing csv to {hp.data.targeted_final.path}")
df_unknown_samples_ids.write_csv(hp.data.targeted_final.path)
# %%
project.log.info("END: preprocessing barcode treated data.")
# %%
# ---
