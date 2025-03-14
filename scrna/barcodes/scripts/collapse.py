
import altair as alt
import helpers as hp
import polars as pl

project = hp.Project("scrna/barcodes", "scrna.barcodes.collapse")

hp.setup_altair()


def bar_with_percent(df, y, y_title):
    base = alt.Chart(df).encode(
        alt.X("length:O").axis(ticks=False, labelAngle=0), alt.Y(y).title(y_title)
    )
    return (
        base.mark_bar() + base.mark_text(yOffset=-5).encode(text=alt.Text("percent"))
    ).properties(width=400)


@project.save_plot()
def unique_barcode_length(df):
    data = (
        df.select("barcode", "length")
        .unique()["length"]
        .value_counts()
        .with_columns(
            (pl.col("count") / pl.col("count").sum() * 100).round(2).alias("percent")
        )
    )
    return bar_with_percent(data, "count", "unique barcodes")


@project.save_plot()
def total_reads_barcode_length(df):
    data = (
        df["length"]
        .value_counts()
        .with_columns(
            (pl.col("count") / pl.col("count").sum() * 100).round(2).alias("percent")
        )
    )
    return bar_with_percent(data, "count", "total reads")


def collapse(barcodes, barcode, length):
    if length > 18:
        return barcode

    for ref in barcodes.filter(pl.col("length") > length)["barcode"].to_list():
        if ref.startswith(barcode) or ref.endswith(barcode):
            return ref

    return barcode


@project.save_plot()
def collapsed_barcode_heatmap(collapsed):
    return (
        alt.Chart(
            collapsed.filter(pl.col("length") != pl.col("collapsed-len"))
            .select("length", "collapsed-len", "barcode", "collapsed-barcode")
            .unique()
            .group_by("length", "collapsed-len")
            .agg(pl.len())
        )
        .mark_rect()
        .encode(
            alt.X("length:O").axis(labelAngle=0).scale(domain=list(range(10, 18))),
            alt.Y("collapsed-len:O")
            .scale(domain=list(range(10, 21)[::-1]))
            .title("collapsed length"),
            color=alt.Color("len").scale(scheme="blues").title(("unique", "barcodes")),
        )
    ).configure_axis(ticks=False)


def main():
    project.log.info("BEGIN: collapsing barcodes.")
    raw = pl.read_csv(project.paths.data / "combined-counts.csv")

    # first we will deduplicate the data to remove PCR bias
    project.log.info("PCR Duplication Bias:")
    project.log.info(
        raw.group_by("umi")
        .agg(pl.col("info").n_unique())["info"]
        .value_counts()
        .sort("count")
    )

    df = raw.select("sample", "cell", "umi", "barcode").unique()
    df = df.with_columns(pl.col("barcode").str.len_chars().alias("length"))

    unique_barcode_length(df)
    total_reads_barcode_length(df)

    # barcodes sorted by length and number of umi's
    # to give collapse preference to prevalent barcodes
    barcodes = (
        df.group_by("barcode", "length")
        .agg(pl.col("umi").n_unique())
        .sort("umi", "length", descending=True)
    )

    collapsed_barcodes = pl.DataFrame(
        (
            (barcode, collapse(barcodes, barcode, length))
            for barcode, length in barcodes.select("barcode", "length").iter_rows()
        ),
        schema=("barcode", "collapsed-barcode"),
    )

    collapsed = df.join(collapsed_barcodes, on="barcode").with_columns(
        pl.col("collapsed-barcode").str.len_chars().alias("collapsed-len")
    )

    collapsed_barcode_heatmap(collapsed)

    project.log.info("Unique umi's and collapsed barcodes")
    project.log.info(
        collapsed.filter(pl.col("length") != pl.col("collapsed-len")).select(
            pl.col("umi").n_unique(), pl.col("collapsed-barcode").n_unique()
        )
    )

    project.log.info(f"writing csv to {project.paths.data / 'collapsed-barcodes.csv'}")
    (
        collapsed.select(
            "sample", "cell", "umi", "barcode", "collapsed-barcode"
        ).write_csv(project.paths.data / "collapsed-barcodes.csv")
    )
    project.log.info("END: collapsing barcodes.")


if __name__ == "__main__":
    main()
