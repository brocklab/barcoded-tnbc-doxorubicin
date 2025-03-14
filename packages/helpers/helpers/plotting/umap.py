import logging

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

# TODO: control how NA is handled
# TODO: add single value plotter with cell_idx + observation

# where the rest of the cells are set to NA
logger = logging.getLogger("helpers.plotting")


class Umap:
    def __init__(self, adata: AnnData, key: str, subset_idx=None):
        self._adata = adata.copy()  # copy?
        if subset_idx is not None:
            self._adata.obs.loc[~self._adata.obs.index.isin(subset_idx), key] = float(
                "nan"
            )
        self.key = key
        self.df = sc.get.obs_df(
            self._adata,
            keys=[key],
            obsm_keys=[("X_umap", 0), ("X_umap", 1)],
        )

    def _add_axis(self, ax):
        arrow_patch_kwargs = dict(
            arrowstyle="-|>,head_width=.15",
            facecolor="black",
            mutation_scale=20,
            transform=ax.transAxes,
            clip_on=False,
        )

        arrY = mpatches.FancyArrowPatch((0, 0), (0, 0.2), **arrow_patch_kwargs)
        arrX = mpatches.FancyArrowPatch((0, 0), (0.2, 0), **arrow_patch_kwargs)
        ax.add_patch(arrY)
        ax.add_patch(arrX)
        ax.annotate(
            "UMAP 2",
            (0, 0.5),
            xycoords=arrY,
            ha="right",
            va="center",
            rotation=90,
            fontsize=9,
        )
        ax.annotate(
            "UMAP 1", (0.5, 0), xycoords=arrX, ha="center", va="top", fontsize=9
        )

    def plot(
        self,
        axis: bool = True,
        save=None,
        cmap=None,
        title=None,
        plotna=True,
        legend=True,
        legend_kwargs=dict(),
        use_uns_colors=False,
    ):
        with plt.rc_context({"legend.markerscale": 8, "legend.frameon": False}):
            fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)

            scatter_kwargs = dict(s=1.5, alpha=0.7, edgecolor="none")

            if self.df[self.key].dtype.name == "category":
                if use_uns_colors and not cmap:
                    cmap = dict(
                        zip(
                            self.df[self.key].cat.categories,
                            self._adata.uns[f"{self.key}_colors"],
                        )
                    )
                if plotna:
                    df_na = self.df.loc[self.df[self.key].isna()]
                    ax.scatter(
                        df_na["X_umap-0"],
                        df_na["X_umap-1"],
                        **scatter_kwargs,
                        c="lightgrey",
                    )

                def plot_cat(df, ax, cat, scatter_kwargs):
                    # pandas seems to be passing empty
                    # dataframes to this function, so let's just bail out
                    if len(df) == 0:
                        return

                    value = df[cat].iloc[0]

                    if cmap:
                        scatter_kwargs.update(dict(color=cmap[value]))

                    ax.scatter(
                        df["X_umap-0"], df["X_umap-1"], **scatter_kwargs, label=value
                    )

                self.df.groupby(self.key, observed=True).apply(
                    plot_cat,
                    ax,
                    self.key,
                    scatter_kwargs,
                )

                if legend:
                    if legend_kwargs:
                        fig.legend(
                            # bbox_to_anchor=(1.04, 0.5), loc="center left",
                            borderaxespad=0,
                            **legend_kwargs,
                        )
                    else:
                        ax.legend()

            else:
                pts = ax.scatter(
                    self.df["X_umap-0"],
                    self.df["X_umap-1"],
                    s=0.7,
                    c=self.df[self.key],
                )
                # Adding the colorbar
                cbaxes = fig.add_axes(
                    # This is the position for the colorbar
                    [1, 0.1, 0.03, 0.8]
                )
                if legend:
                    fig.colorbar(pts, cax=cbaxes)

            if title:
                ax.set_title(title)
            elif title is None:
                ax.set_title(self.key)

            ax.axis("off")

            if axis:
                self._add_axis(ax)

            if save:
                logger.info(f"saving umap to: {save}")
                fig.savefig(save, bbox_inches="tight", dpi=300)

            plt.show()
            plt.close()
