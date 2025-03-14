# TODO: move pycashier rules to their own file
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, List

import polars as pl
import helpers as hp

pro = hp.Project("targeted")
dataset = hp.data.get_sc_dataset("231-1KB3")
plots = [
    pro.paths.outs / f
    for f in [
        "barcodes-post-treatment-per-sample-per-group-swarmplot.svg",
        # "top-surviving-barcodes-dotplot.svg",
        # "top-surviving-barcodes-by-dose.svg",
        "percent-unknown.svg",
        "post-treatment-barcode-dotplot.svg",
        "clone-cluster-unique-barcodes-dots.svg",
        "clone-cluster-abundance-passage-area.svg",
    ]
]


def run_notebook(notebook: str) -> str:
    return f"python targeted/notebooks/{notebook}.py"


@dataclass
class Job:
    id: str
    output: Path
    sample_str: str

    @property
    def extract_outs(self) -> Path:
        return f"data/pycashier-outs-{self.id}"

    @property
    def paired_end_fastqs(self) -> List[Path]:
        return (
            pro.paths.data
            / "raw-paired-end"
            / f"{self.id}-{sample}.R{read}.raw.fastq.gz"
            for read in ("1", "2")
            for sample in self.samples
        )

    @property
    def raw_fastqs(self) -> List[Path]:
        return (
            pro.paths.data / "raw" / f"{self.id}-{sample}.merged.raw.fastq"
            for sample in self.samples
        )

    @property
    def samples(self) -> str:
        return self.sample_str.split(",")

    @property
    def pycashier_samples(self) -> str:
        return ",".join(
            [f"{self.id}-{sample}" for sample in self.sample_str.split(",")]
        )

    @property
    def log_path(self) -> Path:
        return f"logs/pycashier-{self.id}.log"


jobs = [
    Job("JA21116", pro.paths.data / "counts-library.tsv", "libB3"),
    Job(
        "JA21458",
        pro.paths.data / "counts-treated.tsv",
        "1KB3,A1,A2,A3,A4,A5,A6,B1,B2,B4,B5,B6,C2,C3,C4,C5,C6",
    ),
    Job(
        "JA22214",
        pro.paths.data / "counts-outgrowth.tsv",
        "1KB3-P17,1KB3-P19,1KB3-P21,1KB3-P23,1KB3-P25,1KB3-P27,1KB3-P29,1KB3-P31",
    ),
    Job(
        "JA23388",
        pro.paths.data / "counts-esam.tsv",
        "1KB3-0712,1KB3-0719,1KB3-0729,1KB3-810,1KB3-EN-0729,"
        "1KB3-EP-0729,1KB3-EN-0801,1KB3-EN2-0809,1KB3-EP-0801,"
        "1KB3-EP2-0809",
    ),
]


rule pycashier_outs:
    input:
        [job.output for job in jobs],


for job in jobs:

    rule:
        name:
            f"targeted_pycashier_merge_{job.id}"
        input:
            job.paired_end_fastqs,
        output:
            job.raw_fastqs,
        params:
            samples=job.pycashier_samples,
            log_path=job.log_path,
        shell:
            """
          cd targeted
          mkdir -p logs
          pycashier merge \
              --samples {params.samples} \
              --log-file {params.log_path} \
              --yes
          """

    rule:
        name:
            f"targeted_pycashier_{job.id}"
        input:
            job.raw_fastqs,
        output:
            job.output,
        params:
            jobid=job.id,
            samples=job.pycashier_samples,
            log_path=job.log_path,
            extract_outs=job.extract_outs,
        shell:
            """
            cd targeted
            mkdir -p logs
            mkdir -p pycashier/outs/{params.jobid}
            pycashier extract \
                --samples {params.samples} \
                --log-file {params.log_path} \
                --yes \
                --output {params.extract_outs}
            pycashier receipt \
                --samples {params.samples} \
                --input {params.extract_outs} \
                --output {output}
            """


rule preprocess_treated:
    input:
        pro.paths.data / "counts-treated.tsv",
        pro.paths.data / "barcode-bgID.csv",
        pro.paths.data / "DM2103-bfp.tsv",
    output:
        pro.paths.outs / "proportion-barcodes-in-pretreatment.svg",
        pro.paths.data / "treated-filtered.csv",
        hp.data.targeted_final.path,
    shell:
        "python targeted/scripts/preprocess_treated.py"


rule barcode_ids:
    input:
        pro.paths.data / "counts-library.tsv",
    output:
        pro.paths.data / "barcode-bgID.csv",
    shell:
        "python targeted/scripts/barcode-to-bgID.py"


rule supplement:
    input:
        pro.paths.data / "DM2103-bfp.tsv",  # <- raw data file
    output:
        pro.paths.outs / "bfp-positive-sample.svg",
    shell:
        run_notebook("supplement")


rule allowlist:
    input:
        library=pro.paths.data / "counts-library.tsv",
    output:
        allowlist=pro.paths.data / "allowlist.csv",
    run:
        (
            hp.data.load_pycashier_outs(input.library)
            .select("barcode")
            .write_csv(output.allowlist)
        )


rule preprocess_outgrowth:
    input:
        pro.paths.data / "counts-outgrowth.tsv",
        pro.paths.data / "counts-library.tsv",
    output:
        pro.paths.data / "outgrowth.csv",
        pro.paths.outs / "outgrowth-barcodes-in-library.svg",
    shell:
        "python targeted/scripts/preprocess_outgrowth.py"


rule outgrowth:
    input:
        pro.paths.data / "outgrowth.csv",
    output:
        pro.paths.outs / "top-ten-P31-barcodes-abundance-change.svg",
        pro.paths.outs / "unique-barcodes-outgrowth.svg",
    shell:
        run_notebook("outgrowth")


rule subpopulations:
    input:
        pro.paths.data / "counts-esam.tsv",
        pro.paths.data / "barcode-bgID.csv",
        hp.data.obs_final.path,
    output:
        *[
            pro.paths.outs / f
            for f in (
                "esam-separated-totals.svg",
                "clone-cluster-per-sample-subpopulations.svg",
                "frequency-clone-cluster-sorted-bar-dot.svg",
                "subpopulation-dotplot-barcodes.svg",
                "subpopulation-dotplot-top-barcodes.svg",
            )
        ],
    shell:
        run_notebook("subpopulations")


final_figures = [
    pro.paths.outs / f
    for f in [
        "clone-cluster-unique-barcodes-dots.svg",
        "clone-cluster-abundance-passage-area.svg",
        "percent-unknown.svg",
        "barcodes-post-treatment-per-sample-per-group-swarmplot.svg",
        "post-treatment-barcode-dotplot.svg",
        "top-14-barcode-dot-plot-facet.svg",
        "outgrowth-persistence-vs-initial-abundance.svg",
        "shared-bar-post.svg",
        "barcode-dotplot.svg",
        "clone-cluster-per-treated-samples.svg",
    ]
]


rule figures:
    input:
        hp.data.targeted_final.path,
        hp.data.obs_final.path,
        pro.paths.data / "outgrowth.csv",
    output:
        *final_figures,
    shell:
        run_notebook("figures")


rule growth:
    input:
        hp.data.targeted_final.path,
        hp.data.obs_final.path,
        pro.paths.data / "outgrowth.csv",
    output:
        *[       pro.paths.outs / f for f in (

'growth-rates-vs-treated-abundance.svg',
'boxplot-growth-rate-clone-cluster.svg',
'barcode-dotplot-outgrowth.svg',
        )
 ]
    shell:
        run_notebook("growth")
