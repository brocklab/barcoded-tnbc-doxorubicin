import altair as alt
import helpers as hp
import polars as pl
import scanpy as sc

hp.setup_altair()
project = hp.Project("scrna/barcodes", "scrna.barcodes.assign")


def targeted_allowlist():
    # TODO: HARDCODED
    return (
        pl.read_csv(project.paths.root / "targeted" / "data" / "allowlist.csv")[
            "barcode"
        ]
        .unique()
        .to_list()
    )


def y_jitter():
    return (
        alt.Chart()
        .encode(
            alt.Y("in-library").title("found in library").axis(ticks=False, offset=10),
            color=alt.Color("in-library").legend(None),
            yOffset="jitter:Q",
        )
        .transform_calculate(
            # Generate Gaussian jitter with a Box-Muller transform
            jitter="sqrt(-2*log(random()))*cos(2*PI*random())"
        )
    ).properties(
        width=400,
        height=200,
    )


@project.save_plot()
def barcodes_in_library_by_umi_total(barcodes):
    return (
        y_jitter()
        .mark_circle(opacity=0.6)
        .encode(
            alt.X("n_umis").scale(type="log").title("average umis"),
        )
        .properties(
            data=(
                barcodes.group_by("barcode", "in-library").agg(pl.col("n_umis").mean())
            ),
        )
    )


@project.save_plot()
def barcodes_in_library_by_cell(barcodes):
    return (
        y_jitter()
        .mark_circle(opacity=0.6)
        .encode(
            alt.X("cell").scale(type="log").title("cells (N)"),
        )
        .properties(
            data=barcodes.group_by("barcode", "in-library").agg(
                pl.col("cell").n_unique()
            )
        )
    )


@project.save_plot()
def barcodes_in_library_by_cell_and_umi_total(barcodes):
    source = barcodes.group_by("barcode", "in-library").agg(
        pl.col("cell").n_unique(), pl.col("n_umis").mean()
    )
    return (
        alt.Chart(source)
        .mark_circle(opacity=0.6)
        .encode(
            alt.Y("cell").scale(type="log").title("cells (N)"),
            alt.X("n_umis").scale(type="log").title("average umis"),
            alt.Color("in-library").title("in library"),
        )
        .configure_axis(offset=5)
        .properties(
            width=200,
            height=200,
        )
    )


@project.save_plot()
def barcodes_per_cell_heatmap(barcodes):
    base = alt.Chart(
        barcodes.group_by("sample", "cell")
        .agg(pl.col("barcode").n_unique().alias("n_barcodes"))
        .group_by("sample", "n_barcodes")
        .agg(pl.len())
        .with_columns(
            (pl.col("len") / pl.col("len").sum().over("sample") * 100).round(2)
        )
    ).encode(
        alt.X("sample:O").axis(labelAngle=0, orient="top"),
        alt.Y("n_barcodes:O").title("barcodes (N)"),
    )

    hmap = base.mark_rect(binSpacing=0, clip=True).encode(
        color=alt.Size("len")
        .title("percent")
        .legend(gradientLength=150)
        .scale(domain=[0, 100], scheme="yelloworangered"),
    )

    text = base.mark_text().encode(
        text=alt.Text("len"),
        color=alt.condition(alt.datum.len < 50, alt.value("black"), alt.value("white")),
    )

    return (
        (hmap + text)
        .properties(width=400, height=200)
        .configure_view(strokeWidth=1, stroke="black")
        .configure_axis(tickBand="extent", ticks=False, labelPadding=5)
    )


@project.save_plot()
def barcodes_per_cell(barcodes):
    return barcodes_per_cell_heatmap(barcodes)


@project.save_plot()
def barcodes_per_cell_w_cutoff(barcodes, cutoff=3):
    return barcodes_per_cell_heatmap(
        barcodes.filter(pl.col("n_umis") > cutoff)
    ).properties(title=f"umi > {cutoff}")


@project.save_plot()
def doublet_vs_n_cells(multiplets):
    data = (
        multiplets.select("cell", "combo", "doublet-score", "doublet")
        .unique()
        .group_by("combo", "doublet")
        .agg(pl.col("doublet-score").mean(), pl.col("cell").n_unique().alias("n_cells"))
    )

    return (
        alt.Chart(data)
        .mark_circle(size=50, opacity=0.5)
        .encode(
            alt.X("doublet-score").axis(offset=10),
            alt.Y("n_cells")
            .scale(
                type="log",
            )
            .axis(offset=10)
            .title("cells (N) "),
            color=alt.Color("doublet"),
        )
        .properties(width=300)
    )


@project.save_plot()
def n_barcodes_by_doublet_status(multiplets):
    return (
        alt.Chart(
            multiplets.select("sample", "cell", "combo", "doublet", "n_barcodes")
            .unique()
            .group_by("n_barcodes", "doublet", "sample")
            .agg(pl.col("cell").n_unique().alias("n_cells"))
        )
        .mark_circle()
        .encode(
            alt.X("n_barcodes:O").axis(offset=5),
            alt.Y("sample").axis(offset=5),
            size=alt.Size("n_cells"),
            column=alt.Column("doublet"),
        )
        .properties(height=100)
    )


@project.save_plot()
def top_combos_with_doublets(multiplets):
    data = (
        multiplets.filter(
            pl.col("combo").is_in(
                multiplets.group_by("combo")
                .agg(pl.col("cell").n_unique())
                .filter(pl.col("cell") > 5)["combo"]
            )
        )
        .group_by("combo", "doublet")
        .agg(pl.len().alias("n_cells"))
    )
    return (
        alt.Chart(data)
        .mark_bar()
        .encode(
            alt.Y("combo")
            .sort("-x")
            .axis(
                labelLimit=500,
                labelFontSize=8,
                labelFont="monospace",
                ticks=False,
                domain=False,
            )
            .title(None),
            alt.X("n_cells").title("cells (N)").scale(domain=(0, 200)),
            color=alt.Color("doublet").legend(orient="bottom-right"),
        )
    ).properties(height=250)


def get_multiple_integrations(barcodes):
    # multiplet barcodes that exist on their own...
    # NOTE: for now pretend "combo" column doesn't exist

    # probably a cleaner way to do this...
    combos = (
        barcodes.group_by("cell")
        .agg(pl.col("barcode").sort())
        .with_columns(
            pl.col("barcode").list.join("-").alias("combo"),
            #   pl.col('barcode').list.len().alias('combo-size')
        )
        .explode("barcode")
        .group_by("barcode")
        .agg(pl.col("combo").unique().alias("unique-combos"))
        .with_columns(pl.col("unique-combos").list.len().alias("n_combos"))
        .explode("unique-combos")
    )

    return (
        combos.filter(pl.col("n_combos") == 1)
        .group_by("unique-combos")
        .agg(pl.col("barcode"))
        .filter(pl.col("barcode").list.len() > 1)["unique-combos"]
        .to_list()
    )


def get_likely_barcode_doublets(barcodes):
    return (
        barcodes.filter(pl.col("n_barcodes") > 1)
        .group_by("combo")
        .agg(pl.col("doublet").mode())
        .explode("doublet")
        .filter(pl.col("doublet"))["combo"]
        .to_list()
    )


def get_cell_assignments(barcodes, obs):
    return (
        barcodes.group_by("cell", "sample")
        .agg(
            pl.col("barcode").sort(),
            pl.col("n_umis").sum().alias("barcode-umis"),
        )
        .with_columns(
            pl.col("barcode").list.len().alias("n_barcodes"),
            (combo := pl.col("barcode").list.join("-")),
            combo.is_in(get_multiple_integrations(barcodes)).alias(
                "multiple-integration"
            ),
            combo.is_in(get_likely_barcode_doublets(barcodes)).alias("barcode-doublet"),
        )
    ).join(
        obs.select(
            "cell",
            (pl.col("scDblFinder_class") == "doublet").alias("cell-doublet"),
            # pl.col("scDblFinder_score").alias("doublet-score"),
        ),
        on="cell",
        how="left",
    )


@project.save_plot()
def n_barcodes_unique_and_total_cells(cells):
    width = 100
    data = (
        cells.filter(
            ~(
                pl.col("barcode-doublet")
                | pl.col("cell-doublet")
                | pl.col("multiple-integration")
            )
        )
        .group_by("sample", "n_barcodes")
        .agg(
            pl.col("barcode").n_unique().alias("unique-barcodes"),
            pl.col("cell").n_unique().alias("n_cells"),
        )
    )
    base = alt.Chart(data).encode(color=alt.Color("n_barcodes:O").title("barcodes (N)"))
    unique_barcodes = (
        base.mark_bar()
        .encode(
            alt.Y("sample").axis(offset=5),
            alt.X("unique-barcodes").title("unique barcodes"),
        )
        .properties(width=width)
    )

    total_cells = (
        base.mark_bar()
        .encode(
            alt.X("n_cells").title("cells (N)"),
            alt.Y("sample").axis(None),
        )
        .properties(width=width)
    )
    return (
        (unique_barcodes | total_cells)
        .configure_axisY(domain=False, ticks=False)
        .configure_legend(titleFontSize=10)
    )


def main():
    project.log.info("BEGIN: assigning barcodes to single cells.")
    barcodes = (
        (pl.read_csv(project.paths.data / "collapsed-filtered-clustered-barcodes.csv"))
        .group_by("sample", "cell", "barcode")
        .agg(pl.col("umi").n_unique().alias("n_umis"))
    )
    barcodes = barcodes.with_columns(
        pl.col("barcode").is_in(targeted_allowlist()).alias("in-library")
    )
    adata = sc.read(hp.data.get_sc_dataset("231-1KB3").postqc)
    obs = hp.get_obs(adata)
    barcodes = barcodes.join(
        obs.select(
            "cell",
            (pl.col("scDblFinder_class") == "doublet").alias("doublet"),
            pl.col("scDblFinder_score").alias("doublet-score"),
        ),
        on="cell",
    )

    project.log.info("number of total umi's of barcodes at a given length")
    with pl.Config(tbl_rows=50):
        project.log.info(
            barcodes.with_columns(pl.col("barcode").str.len_chars().alias("len"))
            .group_by("len")
            .agg(
                pl.col("n_umis").sum().alias("total umi's"),
                pl.col("barcode").n_unique().alias("unique barcodes"),
            )
            .sort("len")
        )

    barcodes = barcodes.filter(pl.col("barcode").str.len_chars().alias("len") >= 19)
    project.log.info(
        "total number of barcodes with only 1 umi: "
        + str(
            barcodes.group_by("barcode")
            .agg(pl.col("n_umis").sum())
            .filter(pl.col("n_umis") == 1)["barcode"]
            .n_unique(),
        )
    )

    project.log.info("total number of cells with only 1 umi: ")
    project.log.info(
        barcodes.group_by("sample", "cell")
        .agg(pl.col("n_umis").sum())
        .filter(pl.col("n_umis") == 1)
        .group_by("sample")
        .agg(pl.col("cell").n_unique())
    )

    barcodes_in_library_by_cell_and_umi_total(barcodes)
    barcodes_in_library_by_umi_total(barcodes)
    barcodes_in_library_by_cell(barcodes)
    barcodes = barcodes.filter(pl.col("in-library")).drop("in-library")

    barcodes_per_cell(barcodes)
    barcodes_per_cell_w_cutoff(barcodes)

    project.log.info("total cells: ")
    project.log.info(barcodes.group_by("sample").agg(pl.col("cell").n_unique()))

    project.log.info("total cells with umi > 3:")
    project.log.info(
        barcodes.filter(pl.col("n_umis") > 3)
        .group_by("sample")
        .agg(
            (unique := pl.col("cell").n_unique()).alias("unique cells"),
            (unique / pl.col("cell").len() * 100).alias("percent total"),
        )
    )
    n_barcodes_col = pl.col("barcode").n_unique().over("cell").alias("n_barcodes")
    total_umis_col = pl.col("n_umis").sum().over("cell").alias("total_umis")
    barcodes = (
        # filter lowly abundant barcodes as possible noise
        barcodes.filter(pl.col("n_umis") > 3)
        .with_columns(n_barcodes_col, total_umis_col)
        .with_columns(
            (pl.col("total_umis") / pl.col("n_barcodes")).alias("expected_n_umis"),
        )
        # filter the slightly more abundnant barcodes that are also likely to be noise
        .filter(~(pl.col("n_umis") < (pl.col("expected_n_umis") / 2)))
        .with_columns(n_barcodes_col, total_umis_col)
        .drop("expected_n_umis")
    )

    # compute canonical barcode combos
    barcodes = barcodes.join(
        barcodes.group_by("cell").agg(pl.col("barcode").sort().alias("combo")),
        on="cell",
    ).with_columns(pl.col("combo").list.join("-"))

    project.log.info("unique barcodes and barcode combos:")
    project.log.info(barcodes.select(pl.col("barcode", "combo").n_unique()))

    multiplets = barcodes.filter(pl.col("n_barcodes") > 1)
    # singletons = barcodes.filter(pl.col("n_barcodes") == 1)
    project.log.info("cells with more than 1 barcode:")
    project.log.info(multiplets.group_by("sample").agg(pl.col("cell").n_unique()))

    multiplets.select("cell", "n_barcodes").unique()["n_barcodes"].value_counts()

    doublet_vs_n_cells(multiplets)
    n_barcodes_by_doublet_status(multiplets)

    project.log.info(
        'Total number of barcode "combos" in more than 5 cells: '
        + str(
            multiplets.group_by("combo")
            .agg(pl.col("cell").n_unique())
            .filter(pl.col("cell") > 5)
            .shape[0]
        )
    )
    top_combos_with_doublets(multiplets)

    cells = get_cell_assignments(barcodes, obs)
    n_barcodes_unique_and_total_cells(cells)

    project.log.info(f"writing csv to {hp.data.cell_barcodes_final.path}")
    cells.write_csv(hp.data.cell_barcodes_final.path)
    project.log.info("END: assigning barcodes to single cells.")


if __name__ == "__main__":
    main()
