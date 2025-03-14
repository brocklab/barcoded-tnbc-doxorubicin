#!/usr/bin/env python3

import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

DATADIR = Path(__file__).parent.parent / "data"


@dataclass
class Job:
    id: str
    pe: bool
    samples: str

    @property
    def dir(self) -> Path:
        return DATADIR / "jobs" / self.id


JOBS = [
    Job("JA21116", True, "libB3"),
    Job(
        "JA21458",
        True,
        "1KB3,A1,A2,A3,A4,A5,A6,B1,B2,B4,B5,B6,C2,C3,C4,C5,C6",
    ),
    Job(
        "JA22214",
        True,
        "1KB3-P17,1KB3-P19,1KB3-P21,1KB3-P23,1KB3-P25,1KB3-P27,1KB3-P29,1KB3-P31",
    ),
    Job(
        "JA23388",
        True,
        "1KB3-0712,1KB3-0719,1KB3-0729,1KB3-810,1KB3-EN-0729,"
        "1KB3-EP-0729,1KB3-EN-0801,1KB3-EN2-0809,1KB3-EP-0801,"
        "1KB3-EP2-0809",
    ),
]


def seqdat(job: Job):
    cmd = f"seqdat move {job.id} -o {job.dir}"
    if job.pe:
        cmd += " --paired-end"
    cmd += " --samples " + job.samples
    p = subprocess.run(shlex.split(cmd))
    if p.returncode != 0:
        print(f"seqdat failed for job: {job.id} with code: {p.returncode}")
        sys.exit(p.returncode)


def link(job: Job):
    (out_dir := DATADIR / ("raw-paired-end" if job.pe else "raw")).mkdir(
        parents=True, exist_ok=True
    )
    for file in job.dir.iterdir():
        if not (file_link := out_dir / f"{job.id}-{file.name}").is_symlink():
            # pycashier with docker needs a project relative symlink
            (file_link).symlink_to(".." / file.relative_to(DATADIR))


def main():
    (DATADIR / "jobs").mkdir(exist_ok=True)

    for job in JOBS:
        if not job.dir.is_dir():
            print(f"Fetching data for: {job.id} with seqdat")
            seqdat(job)

        print(f"Symlinking data for: {job.id}")
        link(job)


if __name__ == "__main__":
    main()
