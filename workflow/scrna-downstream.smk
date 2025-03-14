from pathlib import Path
from typing import Annotated, List

import helpers as hp
from helpers.data import SingleCellDataSet

dataset = hp.data.get_sc_dataset("231-1KB3")
project = hp.Project("scrna/downstream")


def hp_scrna_sh(module: str, args: str = "") -> str:
    return f"python -m helpers.scrna.{module} {dataset.name} {args}"


def run_notebook(module: str) -> str:
    return f"python scrna/downstream/notebooks/{module}.py"


def outs_paths(*filenames: str, prefix: str | None = None) -> List[Path]:
    outs = (project.paths.outs / prefix) if prefix else project.paths.outs
    # make sure this isn't a generator
    return [outs / f for f in filenames]


def umap_paths(*categories: str) -> List[Path]:
    return outs_paths(*[f"umap-{cat}.png" for cat in categories], prefix="umaps")


umaps = umap_paths(
    "sample",
    "leiden",
    "esam",
    "survivor",
    "group",
    "clone-cluster",
    "clone-cluster-1-pt-dox",
    "clone-cluster-2-pt-dox",
    "survivorship-clone-cluster-1",
    "survivorship-clone-cluster-2",
    "per-sample",
    "top-10-barcodes",
    "sample-pt-dox-30hr-only",
    "survivorship-Pretreatment",
    "total-counts",
    "phase",
    "sample-grps-PT-30HR-250-400",
    "sample-grps-PT-30HR-250",
    "sample-grps-PT-30HR",
    "sample-grps-PT",
)

plots = [
    *umaps,
]


rule final:
    input:
        dataset.final,


rule qc:
    input:
        dataset.raw,
    output:
        dataset.postqc,
        *outs_paths(
            "distplot_counts_231-1KB3.svg",
            "distplot_counts_231-1KB3_post_qc.svg",
            "scatter_counts_231-1KB3.svg",
            "scatter_counts_231-1KB3_post_qc.svg",
            "violin_pct_counts_mt_231-1KB3.svg",
            "violin_pct_counts_mt_231-1KB3_post_qc.svg",
            prefix="qc",
        ),
    shell:
        hp_scrna_sh("qc")


rule normalize:
    input:
        dataset.postqc,
    output:
        dataset.normalized,
        project.paths.outs / "normalize" / "normalized_counts_231-1KB3.svg",
    shell:
        hp_scrna_sh("normalize")


rule cluster:
    input:
        dataset.normalized,
    output:
        dataset.clustered,
        *outs_paths(
            "filter_genes_dispersion_231-1KB3.svg",
            "pca_counts_231-1KB3.svg",
            "umap_counts_pct_mt_231-1KB3.png",
            prefix="cluster",
        ),
    shell:
        hp_scrna_sh("cluster")


rule barcodes:
    input:
        dataset.clustered,
        hp.data.targeted_final.path,
        hp.data.cell_barcodes_final.path,
    output:
        dataset.barcoded,
        *outs_paths(
            "clonal-cluster-per-sample.svg",
            "dendrogram-barcode-cluster.svg",
            prefix="barcodes",
        ),
    shell:
        hp_scrna_sh("barcode", "--threshold 225")


rule final_adata:
    input:
        dataset.barcoded,
        hp.data.targeted_final.path,
    output:
        protected(dataset.final),
        hp.data.obs_final.path,
        # TODO: remove this
        project.paths.data / "pretreatment-esam-expression-by-barcode.csv",
    shell:
        run_notebook("finalize_anndata")


rule umaps:
    input:
        dataset.final,
    output:
        *umaps,
        *outs_paths("total_counts-umap.svg", prefix="umaps"),
    shell:
        run_notebook("umaps")


rule umaps_extra:
    input:
        dataset.final,
    output:
        *umap_paths(
            "bg010-by-sample",
            "bg016-by-sample",
            "bg032-by-sample",
            "top-ten-survivor-clones",
        ),
        *outs_paths(
            "umap_bg009-samples.png",
            "umap_bg010-samples.png",
            "umap_bg016-samples.png",
            "umap_bg029-samples.png",
            "umap_bg032-samples.png",
            "umap_bg046-samples.png",
            "umap_bg053-samples.png",
            "umap_bg065-samples.png",
            "umap_bg074-samples.png",
            "umap_bg149-samples.png",
            prefix="umaps",
        ),
    shell:
        run_notebook("umaps_extra")


rule misc_figures:
    input:
        dataset.final,
    output:
        *outs_paths(
            "correlation_matrix_sample-clone-cluster.svg",
            "post-barcode-totals-by-sample-dotplot.svg",
            prefix="misc",
        ),
    shell:
        run_notebook("misc_figures")


rule umap_gene:
    input:
        dataset.final,
    output:
        *outs_paths("umap-{gene}.png", prefix="umaps"),
    shell:
        "python scrna/downstream/scripts/umap_gene.py --gene {wildcards.gene}"
