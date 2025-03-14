#!/usr/bin/env python3

# fmt: off
# https://swydd.dayl.in/#automagic-snippet
if not((_i:=__import__)("importlib.util").util.find_spec("swydd")or
(_src:=_i("pathlib").Path(__file__).parent/"swydd/__init__.py").is_file()):
  _r=_i("urllib.request").request.urlopen("https://swydd.dayl.in/swydd.py")
  _src.parent.mkdir(exist_ok=True);_src.write_text(_r.read().decode())  # noqa
# fmt: on


import json
from pathlib import Path

from swydd import cli, get, sub, task

ROOT = Path(__file__).parent


def figures():
    return json.loads((ROOT / "workflow/figures.json").read_text())


def startswith_outs(f):
    paths = [
        "targeted/outs/",
        "scrna/downstream/outs/",
        "scrna/barcodes/outs/drug-dosing/outs",
    ]
    for p in paths:
        if f.startswith(p):
            return True

    return False


@task
def check_raw():
    """ensure raw files exist at expected locations"""
    files = (ROOT / "data" / "raw-files.txt").read_text().splitlines()
    missing = [f for f in files if not (ROOT / f).is_file()]
    if missing:
        print("\n".join(missing))


@task
def unknown_outs():
    """check for files not covered by a snakemake rule"""
    files = [f for f in get("snakemake --lu").splitlines() if startswith_outs(f)]
    print("\n".join(files))


@task
def figure_check():
    """check that a rule defines all the final figures"""
    panels = [ROOT / p for _, panels in figures().items() for p in panels]
    sub("snakemake -fn " + " ".join(map(str, panels)))


cli()
