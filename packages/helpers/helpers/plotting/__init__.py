import logging
import sys
from pathlib import Path

import altair as alt

from ..utils import config
from .heatmap import heatmap
from .volcano import volcano

logger = logging.getLogger("helpers.plotting")
logger.setLevel(logging.DEBUG)


# TODO: update for new altair
def setup_altair(config=None, override=dict()):
    alt.themes.register(
        "custom",
        config
        if config
        else lambda: {
            **dict(
                config=dict(
                    view=dict(strokeWidth=0),
                    axis=dict(grid=False, domainColor="black", tickColor="black"),
                ),
            ),
            **override,
        },
    )
    alt.themes.enable("custom")


def get_palette(name, groups=None):
    palette = config["colors"][name]
    if groups:
        palette = {k: v for k, v in palette.items() if k in groups}
    return palette


def get_alt_palette(name, groups=None):
    palette = get_palette(name, groups)
    return dict(domain=list(palette.keys()), range=list(palette.values()))


__ALL__ = ["heatmap", "setup_altair", "palette", "volcano"]
