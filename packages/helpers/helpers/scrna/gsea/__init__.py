import multiprocessing
import sys

import gseapy as gp
import polars as pl

import helpers as hp


def rank_genes(df):
    """given a dataframe with gene, pval and log2fc prerank genes for GSEA"""

    # adust p-vals by sys.float_info.min for ranking to account for p-val ~= 0
    df = df.filter(pl.col("log2fc").is_not_null()).with_columns(
        pl.col("pval") + sys.float_info.min
    )
    return (
        df.with_columns(
            ((pl.col("pval").log10() * -1) * pl.col("log2fc").sign()).alias("score")
        )
        .sort("score", "log2fc", descending=True)
        .select("gene", "score")
        .to_pandas()
        .set_index("gene")
    )


def run_gsea_prerank(
    degs, gene_sets="GO_Biological_Process_2023", prerank_kwargs=dict()
):
    "wrapper for gp.prerank"
    rnk = rank_genes(degs)
    kwargs = dict(
        threads=int(multiprocessing.cpu_count() * 0.25),
        min_size=5,
        max_size=1000,
        permutation_num=10000,
        # permutation_num=1000,
        outdir=None,  # don't write to disk
        seed=231,
        verbose=True,
    )
    kwargs.update(prerank_kwargs)
    return pl.from_pandas(gp.prerank(rnk=rnk, gene_sets=gene_sets, **kwargs).res2d)


def split_gsea_term_column(gsea):
    split = pl.col("Term").str.split(" (")
    return gsea.with_columns(
        pl.col("Term").alias("full-term"),
        split.list.get(0).alias("term"),
        split.list.get(1).str.split(")").list.get(0).alias("GOID"),
    )
