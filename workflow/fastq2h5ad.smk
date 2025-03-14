from pathlib import Path
from typing import Annotated, List

import helpers as hp

# TODO: 1/4 of resources
# configfile: "config.yml"
config["localcores"] = 24
config["localmem"] = 200

dataset = hp.data.get_sc_dataset("231-1KB3")

project = hp.Project("scrna/fastq2h5ad")


# hardcode?
def fastqs(sample):
    return list((project.paths.data / "raw" / sample).iterdir())


def cellranger_outs(sample):
    return [
        project.paths.project / "counts" / sample / "outs" / f
        for f in [
            "possorted_genome_bam.bam",
            "filtered_feature_bc_matrix.h5",
            "filtered_feature_bc_matrix/features.tsv.gz",
            "filtered_feature_bc_matrix/matrix.mtx.gz",
            "filtered_feature_bc_matrix/barcodes.tsv.gz",
        ]
    ]


rule all:
    input:
        dataset.raw,


rule counts:
    input:
        dataset.sample_10X_mtx_paths,


for sample in dataset.sample_ids:

    rule:
        name:
            f"cellranger_{sample}"
        input:
            list((project.paths.data / "raw" / sample).iterdir()),
        output:
            cellranger_outs(sample),
        params:
            sample=sample,
            localcore=24,
            localmem=200,
        shell:
            """
      mkdir -p counts/
      ./cellranger-8.0.1/cellranger count \
          --id={params.sample} \
          --sample={params.sample} \
          --fastqs=./data/raw/{params.sample} \
          --transcriptome=./data/references/refdata-gex-GRCh38-2024-A \
          --localcores={params.localcore} \
          --localmem={params.localmem} \
          --create-bam=true \
          --nosecondary
      """


rule mtx_to_h5ad:
    input:
        dataset.sample_10X_mtx_paths,
    output:
        protected(dataset.raw),
    params:
        name=dataset.name,
    shell:
        "python -m helpers.scrna.loaders.default {params.name}"
