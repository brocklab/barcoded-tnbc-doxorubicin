from pathlib import Path
from typing import Annotated, List

import helpers as hp
from helpers.data import SingleCellDataSet

dataset = hp.data.get_sc_dataset("231-1KB3")
project = hp.Project("scrna/downstream")


def run_notebook(module: str) -> str:
    return f"python scrna/downstream/notebooks/{module}.py"


def run_script(rest: str) -> str:
    return "python scrna/downstream/scripts/" + rest


def outs_paths(*filenames: str, prefix: str | None = None) -> List[Path]:
    outs = (project.paths.outs / prefix) if prefix else project.paths.outs
    # make sure this isn't a generator
    return [outs / f for f in filenames]


rule shuffle_barcodes:
    input:
        dataset.final,
    output:
        project.paths.data / "randomized-barcode-survivorship.csv",
    shell:
        run_notebook("shuffle_barcode_survivorship")


rule pt_surv_diffexp:
    input:
        dataset.final,
    output:
        protected(
            project.paths.data
            / "degs"
            / "deg-{sample}-cluster-{cluster}-survivorship-low-vs-high.csv"
        ),
    shell:
        run_script(
            "pt_surv_diffexp_mixed_effects.py --sample {wildcards.sample} --cluster {wildcards.cluster}"
        )


rule pt_surv_diffexp_random:
    input:
        dataset.final,
        project.paths.data / "randomized-barcode-survivorship.csv",
    output:
        protected(
            project.paths.data
            / "degs"
            / "deg-{sample}-cluster-{cluster}-survivorship-low-vs-high-randomized.csv"
        ),
    shell:
        run_script(
            "pt_surv_diffexp_mixed_effects.py --sample {wildcards.sample} --cluster {wildcards.cluster} --random"
        )


def pt_surv_deg_paths(base, filename):
    return [
        base / filename.format(sample=sample, cluster=cluster, rand=rand)
        for cluster in ("1", "2")
        for sample in ["Pretreatment", "Dox-550-30HR"]
        for rand in ("", "-randomized")
    ]


pt_surv_diffexp_outs = pt_surv_deg_paths(
    project.paths.data / "degs",
    "deg-{sample}-cluster-{cluster}-survivorship-low-vs-high{rand}.csv",
)


rule pt_surv_diffexp_all:
    input:
        *pt_surv_diffexp_outs,


rule pt_surv_diffexp_followup:
    input:
        dataset.final,
        *pt_surv_diffexp_outs,
    output:
        *pt_surv_deg_paths(
            project.paths.outs / "deg-survivorship",
            "pt-survivorship-{sample}-cluster-{cluster}-deg{rand}-volcano.png",
        ),
        *[
            project.paths.outs
            / "deg-survivorship"
            / f"deg-top-genes-violin-{sample}-{cluster}.svg"
            for cluster in ("1", "2")
            for sample in ["Pretreatment", "Dox-550-30HR"]
        ],
    shell:
        run_notebook("pt_surv_diffexp_followup")


rule clone_cluster_deg:
    input:
        dataset.final,
    output:
        protected(project.paths.data / "degs" / "deg-{sample}-clone-cluster-1-vs-2.csv"),
    shell:
        run_script("clone_cluster_diffexp_mixed_effects.py --sample {wildcards.sample}")


clone_cluster_deg_all_outs = [
    project.paths.data / "degs" / f"deg-{sample}-clone-cluster-1-vs-2.csv"
    for sample in ("Pretreatment", "Dox-550-30HR")
]


rule clone_cluster_deg_all:
    input:
        *clone_cluster_deg_all_outs,


rule clone_cluster_deg_followup:
    input:
        dataset.final,
        *clone_cluster_deg_all_outs,
    output:
        outs_paths(
            "heatmap-pretreament-clone-cluster-1-vs-2.svg", prefix="deg-clone-cluster"
        ),
    shell:
        run_notebook("clone_cluster_deg_followup")


rule pt_30hr_deg:
    input:
        dataset.final,
    output:
        protected(
            project.paths.data
            / "degs"
            / "deg-pretreatment-vs-30HR-clone-cluster-{cluster}.csv",
        ),
    shell:
        run_script("pt_30hr_diffexp_mixed_effects.py --cluster {wildcards.cluster}")


pt_30hr_deg_all_outs = [
    project.paths.data
    / "degs"
    / f"deg-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
    for cluster in ("1", "2")
]


rule pt_30hr_deg_all:
    input:
        *pt_30hr_deg_all_outs,


rule pt_30hr_gsea:
    input:
        project.paths.data
        / "degs"
        / "deg-pretreatment-vs-30HR-clone-cluster-{cluster}.csv",
    output:
        protected(
            project.paths.data
            / "gsea"
            / "gsea-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
        ),
    shell:
        run_script("pt_30hr_gsea.py --cluster {wildcards.cluster}")


pt_30hr_gsea_all_outs = [
    project.paths.data
    / "gsea"
    / f"gsea-pretreatment-vs-30HR-clone-cluster-{cluster}.csv"
    for cluster in ("1", "2")
]


rule pt_30hr_gsea_all:
    input:
        *pt_30hr_gsea_all_outs,


rule pt_30hr_deg_followup:
    input:
        dataset.final,
        *pt_30hr_deg_all_outs,
        *pt_30hr_gsea_all_outs,
    output:
        *outs_paths(
            "pt-vs-30hr-clone-clusters-gsea.svg",
            "pt-vs-30hr-clone-cluster-deg-correlation.svg",
            "heatmap-pt-vs-30hr-clone-clusters-deg-up250-down250.svg",
            "pt-vs-30hr-clone-cluster-1-deg-volcano.png",
            "pt-vs-30hr-clone-cluster-2-deg-volcano.png",
            prefix="pt-30hr",
        ),
    shell:
        run_notebook("pt_30hr_diffexp_followup")


rule new_populations_deg:
    input:
        dataset.final,
    output:
        protected(
            project.paths.data
            / "degs"
            / "deg-new-population-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
        ),
    shell:
        run_script(
            "new_populations_diffexp_mixed_effects.py --sample {wildcards.sample} --bgID {wildcards.bgID}"
        )


new_populations = [
    ("Dox-250-1", "bg009"),
    ("Dox-250-1", "bg032"),
    ("Dox-250-4", "bg016"),
    ("Dox-250-6", "bg074"),
    ("Dox-250-6", "bg032"),
    ("Dox-400-2", "bg053"),
    ("Dox-400-4", "bg149"),
    ("Dox-400-5", "bg016"),
]


new_populations_deg_all_outs = [
    project.paths.data
    / "degs"
    / f"deg-new-population-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    for sample, bgID in new_populations
]


rule new_populations_deg_all:
    input:
        *new_populations_deg_all_outs,


rule new_populations_deg_followup:
    input:
        dataset.final,
        *new_populations_deg_all_outs,
    output:
        *outs_paths(
            "bar-newly-upregulated-totals.svg",
            "dotplot_new-upregulated-genes-post-treatment-clones.svg",
            *[
                f"umap-{sample}-{bgID}-Pretreatment-cluster.png"
                for sample, bgID in new_populations
            ],
            prefix="new-populations",
        ),
    shell:
        run_notebook("new_populations_diffexp_followup")


rule divergent_recovery_deg:
    input:
        dataset.final,
    output:
        protected(
            project.paths.data
            / "degs"
            / "deg-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
        ),
    shell:
        "python scrna/downstream/scripts/divergent_recovery_diffexp.py --sample {wildcards.sample} --bgID {wildcards.bgID}"


divergent_bgIDs = [
    "bg032",
    "bg016",
    "bg010",
]

divergent_bgIDs_samples = [
    ("bg032", "Dox-250-1"),
    ("bg032", "Dox-250-6"),
    ("bg016", "Dox-250-4"),
    ("bg016", "Dox-400-5"),
    ("bg010", "Dox-250-4"),
    ("bg010", "Dox-400-4"),
]


divergent_recovery_deg_all_outs = [
    project.paths.data
    / "degs"
    / f"deg-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    for bgID, sample in divergent_bgIDs_samples
]


rule divergent_recovery_deg_all:
    input:
        *divergent_recovery_deg_all_outs,


rule divergent_recover_gsea:
    input:
        project.paths.data
        / "degs"
        / "deg-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv",
    output:
        protected(
            project.paths.data
            / "gsea"
            / "gsea-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
        ),
    shell:
        run_script(
            "divergent_recovery_gsea.py --sample {wildcards.sample} --bgID {wildcards.bgID}"
        )


divergent_recovery_gsea_all_outs = [
    project.paths.data
    / "gsea"
    / f"gsea-divergent-recovery-{sample}-bgID-{bgID}-vs-Pretreatment.csv"
    for bgID, sample in divergent_bgIDs_samples
]


rule divergent_recovery_gsea_all:
    input:
        *divergent_recovery_gsea_all_outs,


rule divergent_recovery_deg_followup:
    input:
        *divergent_recovery_deg_all_outs,
        *divergent_recovery_gsea_all_outs,
    output:
        *outs_paths(
            *[
                f"overlapping-genes-{bgID}-pre-vs-post-treatment.svg"
                for bgID in divergent_bgIDs
            ],
            *[
                f"divergent-barcode-{bgID}-{sample}-deg-volcano.png"
                for bgID, sample in divergent_bgIDs_samples
            ],
            *[f"gsea-top5-significant-{bgID}.svg" for bgID in divergent_bgIDs],
            prefix="divergent",
        ),
    shell:
        run_notebook("divergent_recovery_diffexp_followup")
