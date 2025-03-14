import scanpy as sc


def setup_scanpy(outs):
    sc.settings.figdir = outs
    sc.settings.verbosity = 3
