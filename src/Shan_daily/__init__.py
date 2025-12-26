"""
Daily-scale drought, heatwave, and compound event identification methods (Python implementation).

Main entry points:
- identify_compound_events: Complete workflow from raw precipitation/temperature daily data to SPI/SHI, thresholds, events, and compound results.
- identify_extremes: Single-variable extreme event identification workflow.
"""

from .daily_2_events import daily_2_events
from .compound4types import identify_compound
from .pipeline import identify_compound_events, identify_extremes
from .PRM_extreme_identification import PRM_extreme_identification
from .remo_merg import remo_merg
from .SI_best_distribution import SI_best_distribution
from .SI_nonparametric import SI_nonparametric
from .plot_extremes_identification import plot_extremes_identification
from .to_365 import to_365

__all__ = [
    "identify_compound_events",
    "identify_extremes",
    "identify_compound",
    "PRM_extreme_identification",
    "remo_merg",
    "daily_2_events",
    "SI_nonparametric",
    "SI_best_distribution",
    "plot_extremes_identification",
    "to_365",
]

