import logging
from multiprocessing import cpu_count

import anndata as ad
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

rcb.logger.setLevel(logging.DEBUG)


logger = logging.getLogger("helpers.scrna.qc.doublets")


def run_scDblFinder(adata: ad.AnnData, cores=None) -> None:
    """
    Args:
        adata: anndata with count matrix to determine doublets
        cores: number of cores to use (default: 1/4 of available)
    """
    if not cores:
        cores = int(cpu_count() / 4)

    data_mat = adata.X.T
    samples = adata.obs["sample"]
    logger.debug("loading R packages")
    importr("SingleCellExperiment")
    importr("scDblFinder")
    importr("BiocParallel")

    detect_doublets = ro.functions.wrap_r_function(
        ro.r(
            """
function(data_mat, samples, cores){
    print("data loaded, proceeding with doublet detection")
    bp <- MulticoreParam(cores, RNGseed=9)
    sce = scDblFinder(
        SingleCellExperiment(
            list(counts=data_mat),
        ),
        samples = samples,
        BPPARAM=bp
    )
    return(
        list(
            score=sce$scDblFinder.score,
            class=sce$scDblFinder.class
        )
    )
}
"""
        ),
        "detect_doublets",
    )

    logger.debug("running scDblFinder")
    with (
        ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
    ).context():
        outs = detect_doublets(data_mat, samples, cores)

    adata.obs["scDblFinder_score"] = outs["score"]
    adata.obs["scDblFinder_class"] = outs["class"]
