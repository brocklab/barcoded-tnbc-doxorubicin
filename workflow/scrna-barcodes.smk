import helpers as hp
import polars as pl

project = hp.Project("scrna/barcodes")
config["threads"] = config.get("threads", 20)
dataset = hp.data.get_sc_dataset("231-1KB3")

cellranger_counts = project.paths.root / "scrna/fastq2h5ad/counts"
sam_paths = (
    [
        project.paths.data / "sams" / f"{sample}.unmapped.sam"
        for sample in dataset.sample_ids
    ],
)


rule all:
    input:
        hp.data.cell_barcodes_final.path,


rule sams:
    input:
        sam_paths,


for sample in dataset.sample_ids:

    rule:
        name:
            f"make_sam_{sample}"
        input:
            cellranger_counts / sample / "outs/possorted_genome_bam.bam",
        output:
            project.paths.data / "sams" / f"{sample}.unmapped.sam",
        shell:
            "samtools view -@ {config[threads]} -f 4 {input} > {output}"


pycashier_outs = [
    project.paths.data / "pycashier/outs" / f"{sample}.umi_cell_labeled.barcode.tsv"
    for sample in dataset.sample_ids
]


rule pycashier:
    input:
        *sam_paths,
    output:
        *pycashier_outs,
    shell:
        """
      cd scrna/barcodes
      mkdir -p data/pycashier
      pycashier scrna -y
      """


rule combine:
    input:
        *pycashier_outs,
    output:
        project.paths.data / "combined-counts.csv",
    run:
        pl.concat(
            (
                pl.read_csv(f, separator="\t").with_columns(
                    sample=pl.lit(f.name.split(".")[0])
                )
                for f in map(Path, input)
            )
        ).write_csv(output[0])


def to_outs(*files: List[str]) -> List[Path]:
    return (project.paths.outs / f for f in files)


rule collapse:
    input:
        project.paths.data / "combined-counts.csv",
    output:
        project.paths.data / "collapsed-barcodes.csv",
        *to_outs(
            "unique-barcode-length.svg",
            "total-reads-barcode-length.svg",
            "collapsed-barcode-heatmap.svg",
        ),
    shell:
        "python scrna/barcodes/scripts/collapse.py"


rule filter:
    input:
        project.paths.data / "collapsed-barcodes.csv",
        dataset.postqc,
    output:
        project.paths.data / "collapsed-filtered-barcodes.csv",
        project.paths.data / "barcodes-pre-clustering.tsv",
        *to_outs(
            "filtered-cells-bar.svg",
        ),
    shell:
        "python scrna/barcodes/scripts/filter.py"


rule starcode:
    input:
        project.paths.data / "barcodes-pre-clustering.tsv",
    output:
        project.paths.data / "barcodes-post-clustering.tsv",
    shell:
        "starcode -d 3 -r 3 --print-clusters -i {input} -o {output}"


rule post_cluster:
    input:
        project.paths.data / "barcodes-post-clustering.tsv",
        project.paths.data / "collapsed-filtered-barcodes.csv",
    output:
        project.paths.data / "collapsed-filtered-clustered-barcodes.csv",
    shell:
        "python scrna/barcodes/scripts/post_cluster.py"


rule assign:
    input:
        project.paths.root / "targeted" / "data" / "allowlist.csv",
        project.paths.data / "collapsed-filtered-clustered-barcodes.csv",
        dataset.postqc,
    output:
        hp.data.cell_barcodes_final.path,
        *to_outs(
            "barcodes-in-library-by-umi-total.svg",
            "barcodes-in-library-by-cell.svg",
            "barcodes-in-library-by-cell-and-umi-total.svg",
            "barcodes-per-cell-heatmap.svg",
            "barcodes-per-cell.svg",
            "barcodes-per-cell-w-cutoff.svg",
            "doublet-vs-n-cells.svg",
            "n-barcodes-by-doublet-status.svg",
            "top-combos-with-doublets.svg",
            "n-barcodes-unique-and-total-cells.svg",
        ),
    shell:
        "python scrna/barcodes/scripts/assign.py"
