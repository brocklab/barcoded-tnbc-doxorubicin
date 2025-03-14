#!/usr/bin/env python3

from collections import namedtuple
from pathlib import Path

SAMPLES = ["PT", "Dox", "A4", "A6"]
ROOT = Path(__file__).parent.parent
DATA_PATH = ROOT / "data/raw"

DataSet = namedtuple("DataSet", "src_path samples renames", defaults=[None, [], dict()])

datasets = [
    DataSet(
        Path("/stor/work/Brock/sequencing/db/JA21437/data"),
        ["PT", "Dox", "A4", "A6"],
        {"PT": "PT-1"},
    ),
    DataSet(
        Path("/stor/work/Brock/sequencing/db/JA24256/data"),
        ["PT", "A1", "B2", "B5", "B6"],
        {"PT": "PT-2"},
    ),
]


def link_sample(src_path, sample, final_name):
    (DATA_PATH / final_name).mkdir(exist_ok=True, parents=True)
    for f in src_path.glob(f"*{sample}*/*"):
        target = DATA_PATH / final_name / f.name.replace(sample, final_name)
        print(f, "-->", target)
        target.symlink_to(f)


def main():
    for dataset in datasets:
        for sample in dataset.samples:
            final_name = dataset.renames.get(sample, sample)
            link_sample(dataset.src_path, sample, final_name)


if __name__ == "__main__":
    main()
