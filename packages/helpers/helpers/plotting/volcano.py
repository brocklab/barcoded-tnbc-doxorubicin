import logging
import sys
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

logger = logging.getLogger("helpers.plotting")


def volcano(
    dedf: pl.DataFrame,
    pval_cutoff=0.05,
    logFC_cutoff=0.25,
    save=None,
    save_kwargs: Dict[Any, Any] = dict(),
):
    """
    Arguments:
    dedf: pl.DataFrame containing (gene, pval, log2fc)
    pval_cutoff: where to draw horizontal dash
    logFC_cutoff: where to draw vertical dash
    save: str/Path: save file path
    save_kwargs: Dict[Any,Any]: params for fig.savefig


    Returns:
    matplotlib figure
    """
    if dedf["pval"].min() == 0:
        logger.info("zero p-values exist, adjusting with float minimum")
        dedf = dedf.with_columns(pl.col("pval") + sys.float_info.min)

    fig, ax = plt.subplots()

    line_kwargs = dict(color="black", linestyle="--")

    sig_genes = pl.col("pval") < pval_cutoff
    up_genes = (pl.col("log2fc") > logFC_cutoff) & sig_genes
    down_genes = (pl.col("log2fc") < -1 * logFC_cutoff) & sig_genes
    else_genes = ~(up_genes | down_genes)

    for f, color in ((up_genes, "blue"), (down_genes, "red"), (else_genes, "grey")):
        ax.scatter(
            x=dedf.filter(f)["log2fc"],
            y=dedf.filter(f)["pval"].log10() * -1,
            c=color,
            s=5,
            alpha=0.8,
        )
    ax.annotate(
        f"{dedf.filter(up_genes).shape[0]} up-regulated",
        (0.95, 1.01),
        xycoords="axes fraction",
        ha="right",
    )
    ax.annotate(
        f"{dedf.filter(down_genes).shape[0]} down-regulated",
        (0.05, 1.01),
        xycoords="axes fraction",
        ha="left",
    )

    x_lim = dedf["log2fc"].abs().max() * 1.1
    ax.set_xlim(-1 * x_lim, x_lim)
    ax.axhline(-1 * np.log10(pval_cutoff), **line_kwargs)
    ax.axvline(logFC_cutoff, **line_kwargs)
    ax.axvline(-1 * logFC_cutoff, **line_kwargs)
    ax.set(xlabel=r"$\log_{2}(FC)$", ylabel=r"$-\log_{10}(pval)$")
    if save:
        fig.savefig(save, **save_kwargs)
    return fig
