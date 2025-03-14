import logging
from pathlib import Path

import anndata2ri
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import rpy2.robjects as ro
import scanpy as sc
from rpy2.robjects.packages import importr

import helpers as hp

from ...plotting import heatmap as hmap

importr("MAST")


logger = logging.getLogger("helpers.scrna.diffexp")


# TODO: Deprecate DiffExp class

# class DiffExp:
#     def __init__(self, adata, groupby, group, min_cells=None):
#         self._adata = adata.copy()
#         self.group = group
#         self.groupby = groupby
#         if not min_cells:
#             min_cells = int(self._adata.obs.shape[0] * 0.1)
#         sc.pp.filter_genes(self._adata, min_cells=min_cells)
#         self._run_diff_exp(groupby, group)

#     def _run_diff_exp(self, groupby, group):
#         logger.info(f"running wilcoxon rank_genes_groups test on {group} in {groupby}") sc.tl.rank_genes_groups(
#             self._adata,
#             groupby=groupby,
#             method="wilcoxon",
#             use_raw=False,
#             layer="log1p_norm",
#             key_added="wilcoxon",
#             pts=True,
#         )
#         # extract scores
#         dedf = (
#             dedf := sc.get.rank_genes_groups_df(self._adata, group, key="wilcoxon")
#         ).loc[(dedf["pct_nz_group"] > 0.25)]
#         t_stats = (
#             # Get dataframe of DE results for condition vs. rest
#             dedf  # .loc[
#             # (dedf["pct_nz_group"] > 0.25) & (dedf['pct_nz_reference'] < 0.5)
#             # ]
#             # Subset to highly variable genes
#             # .loc[dedf['names'].isin(hvgs.loc[hvgs["highly_variable"]].index)]
#             .set_index("names")[
#                 # Format for decoupler
#                 ["scores"]
#             ]
#             # .sort_values("scores", key=np.abs, ascending=False)
#             .rename_axis([groupby + "-" + group], axis=1)
#         )
#         self.dedf = pl.from_pandas(dedf)
#         self._t_stats = t_stats

#     # def _run_gsea(
#     #     self,
#     # ):
#     #     gsea_results = dc.get_gsea_df(
#     #         self._t_stats,
#     #         "scores",
#     #         reactome_gsea_genesets,
#     #         source="geneset",
#     #         target="genesymbol",
#     #     )

#     #     self.gsea = pl.from_pandas(gsea_results)
#     #     # return pl.from_pandas(dedf), t_stats, pl.from_pandas(gsea_results)

#     def leading_edge(self, term):
#         if term not in self.gsea["Term"]:
#             raise ValueError(f"{term} not found in gsea results")

#         return (
#             self.gsea.filter(pl.col("Term") == term)
#             .select("Leading edge")
#             .item()
#             .split(";")
#         )

#     def volcano_plot(self, pval_cutoff=0.05, logFC_cutoff=0.25, save=None):
#         dedf = self.dedf
#         fig, ax = plt.subplots()

#         line_kwargs = dict(color="black", linestyle="--", alpha=0.5)

#         sig_genes = pl.col("pvals") < pval_cutoff
#         up_genes = (pl.col("logfoldchanges") > logFC_cutoff) & sig_genes
#         down_genes = (pl.col("logfoldchanges") < -1 * logFC_cutoff) & sig_genes
#         else_genes = ~(up_genes | down_genes)

#         for f, color in ((up_genes, "blue"), (down_genes, "red"), (else_genes, "grey")):
#             ax.scatter(
#                 x=dedf.filter(f)["logfoldchanges"],
#                 y=dedf.filter(f)["pvals"].log10() * -1,
#                 c=color,
#                 s=5,
#                 alpha=0.8,
#             )
#         ax.annotate(
#             f"{dedf.filter(up_genes).shape[0]} up-regulated",
#             (0.95, 1.01),
#             xycoords="axes fraction",
#             ha="right",
#         )
#         ax.annotate(
#             f"{dedf.filter(down_genes).shape[0]} down-regulated",
#             (0.05, 1.01),
#             xycoords="axes fraction",
#             ha="left",
#         )

#         x_lim = dedf["logfoldchanges"].abs().max() * 1.1
#         ax.set_xlim(-1 * x_lim, x_lim)
#         ax.axhline(-1 * np.log10(pval_cutoff), **line_kwargs)
#         ax.axvline(logFC_cutoff, **line_kwargs)
#         ax.axvline(-1 * logFC_cutoff, **line_kwargs)
#         ax.set(xlabel="$\log_{2}(FC)$", ylabel="$\log_{-10}(pval)$")

#         if save:
#             fig.savefig(save)

#         return fig

#     def get_sig_genes(self, logFC_cutoff=0.25, pval_cutoff=0.05, max_genes=None):
#         sig_genes = pl.col("pvals") < pval_cutoff
#         up_genes_filter = (pl.col("logfoldchanges") > logFC_cutoff) & sig_genes
#         down_genes_filter = (pl.col("logfoldchanges") < -1 * logFC_cutoff) & sig_genes
#         # else_genes = ~(up_genes | down_genes)

#         up_genes = (
#             self.dedf.filter(up_genes_filter)
#             .sort("logfoldchanges", descending=True)["names"]
#             .to_list()
#         )

#         down_genes = (
#             self.dedf.filter(down_genes_filter)
#             .sort("logfoldchanges")["names"]
#             .to_list()
#         )

#         return {
#             "up": up_genes[: (max_genes or len(up_genes) + 1)],
#             "down": down_genes[: (max_genes or len(down_genes) + 1)],
#         }

#     def heatmap(
#         self, adata=None, genes=None, n_genes=None, groupby=None, save=None, **kwargs
#     ):
#         if adata is None:
#             adata = self._adata
#         if genes is None:
#             genes = self.get_sig_genes(max_genes=n_genes or 25)
#         if groupby is None:
#             groupby = self.groupby
#         kwargs.update(dict(groupby=groupby))

#         hmap(adata=adata, genes=genes, save=save, **kwargs)


def find_de_MAST_barcode(adata_, obs, reference):
    """
    mixed model with random effects from barcode
    """
    adata_ = adata_.copy()
    categories = adata_.obs[obs].unique().tolist()
    print(f"{obs=}, {categories=}")
    # these asserts should make group impossible to generate incorrectly
    assert len(categories) == 2, "categories an unexpected length"
    assert reference in categories
    group = ""
    try:
        group = "label" + [cat for cat in categories if cat != reference][0]
    except IndexError:
        print(f"expected a second category besides: {reference}")
        print(f"found {len(categories)} categories: {categories}")
    adata_.obs["label"] = adata_.obs[obs]
    sc.pp.filter_genes(adata_, min_cells=adata_.shape[0] * 0.1)
    find_de_MAST_impl = ro.functions.wrap_r_function(
        ro.r((Path(__file__).parent / "mast_mixed_effects_barcode.R").read_text()),
        "fit_model",
    )
    with (
        ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
    ).context():
        degs = find_de_MAST_impl(adata_, group, reference)

    return pl.from_pandas(degs).rename(
        {
            "primerid": "gene",
            "Pr(>Chisq)": "pval",
            "coef": "log2fc",
        }
    )


def find_de_MAST(adata_, obs, reference):
    adata_ = adata_.copy()
    categories = adata_.obs[obs].unique().tolist()
    print(categories)
    # these asserts should make group impossible to generate incorrectly
    assert len(categories) == 2
    assert reference in categories
    try:
        group = "label" + [cat for cat in categories if cat != reference][0]
    except IndexError:
        print(f"expected a second category besides: {reference}")
        print(f"found {len(categories)} categories: {categories}")
    adata_.obs["label"] = adata_.obs[obs]
    sc.pp.filter_genes(adata_, min_cells=adata_.shape[0] * 0.1)
    find_de_MAST_impl = ro.functions.wrap_r_function(
        ro.r((Path(__file__).parent / "mast.R").read_text()),
        "fit_model",
    )
    with (
        ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
    ).context():
        degs = find_de_MAST_impl(adata_, group, reference)

    return pl.from_pandas(degs).rename(
        {
            "primerid": "gene",
            "Pr(>Chisq)": "pval",
            "coef": "log2fc",
        }
    )
