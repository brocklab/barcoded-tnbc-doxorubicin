import inspect
import logging
import shlex
import subprocess as subp
from functools import wraps
from logging import StreamHandler
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Any

import scanpy as sc
from ruamel.yaml import YAML

yaml = YAML(typ="safe")  # default, if not specfied, is 'rt' (round-trip)


ROOT_PATH = Path(__file__).parent.parent.parent.parent.parent


def get_root_dir(d=Path):
    """walk up the directories to find the root directory"""
    # d = Path.cwd()
    if not (d.is_file() or d.is_dir()):
        raise ValueError(f"{d} does not exist.")
    while not ((d / "README.md").is_file() or (d / ".git").is_dir()):
        d = d.parent
    return d


class DirPaths:
    def __repr__(self):
        return (
            f"DirPaths(\n  root={repr(self.root)},\n  project={repr(self.project)},\n)"
        )

    def __init__(
        self,
        start: str | Path | None = None,
        ensure: bool = False,
        outs: str | None = None,
    ) -> None:
        """generate class with project paths"""

        if not (p := ROOT_PATH / start).is_dir():
            raise ValueError(f"expected path at {p} does not exist")

        self._root = ROOT_PATH
        project = get_root_dir(p)
        self._project = project
        self._outs = project / (outs if outs else "outs")
        self._data = project / "data"

        if ensure:
            self.touch()

    def __getattr__(self, name: str):
        if name == "root":
            return self._root
        elif name in ("project", "outs", "data", "logs"):
            if (v := self.__dict__.get(f"_{name}")) is None:
                raise AttributeError(
                    f"paths.{name} not set. Did you use paths.update or set_project?"
                )
            else:
                return v
        else:
            raise AttributeError(f"{name} is not an attribute of DirPaths")

    def __setattr__(self, name: str, value: Any):
        # named_attrs = ("root", "project", "outs", "data", "logs")
        named_attrs = ("root", "project", "outs", "data")
        if name in named_attrs:
            raise AttributeError("Attempted to modify read-only attribute: {name}")
        elif name not in map(lambda x: f"_{x}", named_attrs):
            raise AttributeError(f"{name} is not an attribute of DirPaths")
        else:
            self.__dict__[name] = value

    def touch(self) -> None:
        self.outs.mkdir(exist_ok=True)
        self.data.mkdir(exist_ok=True)

    def update(self, start: Path, ensure: bool = True) -> None:
        project = get_root_dir(start)
        self._project = project
        self._outs = project / "outs"
        self._data = project / "data"

        if ensure:
            self.touch()


class Project:
    def __init__(
        self, path: str, logger: str | None = None, outs: str | None = None
    ) -> None:
        self.paths = DirPaths(path, outs=outs)
        self.paths.touch()
        # TODO: ensure logging handler isn't added multiple times..
        if logger:
            self.log = logging.getLogger(logger)
        else:
            self.log = logging.getLogger(path.replace("/", "."))
        self.log = setup_logger(self.log)

    # TODO: add path kwarg to wrapped function so at call time can be modified
    def save_plot(self, path: Path | None = None, fmt: str = "svg", **save_kwargs):
        def inner(f):
            @wraps(f)
            def wrapped(*args, **kwargs):
                filepath = (
                    path
                    if path
                    else (self.paths.outs / f"{f.__name__.replace('_', '-')}.{fmt}")
                )
                self.log.info(f"saving plot to: {filepath}")
                p = f(*args, **kwargs)
                p.save(filepath, **save_kwargs)
                return p

            return wrapped

        return inner

    def setup_scanpy(self, figdir: Path = None) -> None:
        if not figdir:
            sc.settings.figdir = self.paths.outs
        else:
            sc.settings.figdir = figdir
        sc.settings.verbosity = 3


def setup_logger(logger):
    logger.setLevel(logging.DEBUG)

    if logger.handlers:
        return logger

    formatter = logging.Formatter("%(asctime)s|%(name)s|%(levelname)s|%(message)s")
    sh = StreamHandler()
    if not (ROOT_PATH / "logs").is_dir():
        (ROOT_PATH / "logs").mkdir(exist_ok=True)
    fh = RotatingFileHandler(
        ROOT_PATH / "logs" / "helpers.log", maxBytes=100000, backupCount=10
    )

    for h in (sh, fh):
        h.setFormatter(formatter)
        logger.addHandler(h)

    return logger


def dataset_from_path(path: Path) -> str:
    s = path.name.split("_", 1)
    if len(s) == 1:
        raise ValueError("expected filename to have underscore delimiting dataset name")

    return s[0]


def subprocess(cmd: str, **kwargs: Any) -> subp.CompletedProcess:
    return subp.check_call(shlex.split(cmd), **kwargs)


config = yaml.load(ROOT_PATH / "config.yml")


def reload_config():
    global config
    config = yaml.load(ROOT_PATH / "config.yml")
