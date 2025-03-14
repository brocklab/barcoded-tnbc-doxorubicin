import logging
from datetime import datetime

from ..utils import paths

logger = logging.getLogger("helpers.reporter")
logger.setLevel(logging.DEBUG)


class Reporter:
    def __init__(self, filename: str, title: str = None) -> None:
        self.fd = paths.outs / filename
        self.title = (
            title
            if title
            else self.filename.split(".", 1)[0].replace("-", " ").captalize()
        )

        self._init_report()

    def _init_report(self) -> None:
        if self.fd.is_file():
            logger.info(f"overwriting existing report file: {self.fd}")

        self.fd.unlink()
        self.fd.touch()

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.write(f"""
---
date: {now}
---

# {self.title}
""")

    def write(self, msg: str) -> None:
        with self.fd.open("a") as f:
            f.write(f"{msg}\n")
