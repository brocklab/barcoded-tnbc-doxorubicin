import logging

import matplotlib.pyplot as plt
import scanpy as sc

logger = logging.getLogger("helpers.plotting")


def heatmap(adata, genes, save=None, **kwargs):
    heatmap_kwargs = dict(layer="scaled", vmin=-2, vmax=2, cmap="RdBu_r")
    heatmap_kwargs.update(**kwargs)

    if heatmap_kwargs.get("layer") == "scaled":
        logger.info("calculating scaled layer by running sc.pp.scale on .X")
        adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X

    with plt.rc_context():
        sc.pl.heatmap(adata, var_names=genes, **heatmap_kwargs, show=False)
        if save:
            plt.savefig(save, bbox_inches="tight", dpi=300)
