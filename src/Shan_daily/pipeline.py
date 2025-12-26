from __future__ import annotations

from typing import Any, Literal, Sequence, Tuple

import numpy as np
import pandas as pd

from .PRM_extreme_identification import (
    NEGATIVE_TYPES,
    POSITIVE_TYPES,
    PRM_extreme_identification,
)
from .SI_best_distribution import SI_best_distribution
from .SI_nonparametric import SI_nonparametric
from .daily_2_events import daily_2_events
from .compound4types import identify_compound
from .remo_merg import remo_merg
from .to_365 import to_365

SiMethod = Literal["nonparametric", "best_distribution"]


def _get_thresholds(df: pd.DataFrame) -> Tuple[float, float]:
    """Extract removal/merging thresholds from remo_merg result's first row, return NaN if missing."""
    if df is None or df.empty:
        return float("nan"), float("nan")
    row = df.iloc[0]
    return float(row.get("removal_threshold", np.nan)), float(row.get("merging_threshold", np.nan))


def _sanitize_thresholds(removal: float | None, merging: float | None, default_merging: float) -> tuple[float | None, float | None]:
    """Convert NaN to None/default values for 
    PRM_extreme_identification processing."""
    rem_nan = removal is None or (isinstance(removal, (float, np.floating)) and np.isnan(removal))
    mer_nan = merging is None or (isinstance(merging, (float, np.floating)) and np.isnan(merging))
    rem = None if rem_nan else float(removal)
    mer = default_merging if mer_nan else float(merging)
    return rem, mer


def _default_threshold(extreme_type: str) -> tuple[float, float]:
    """Return default pre-identification thresholds based on 
    extreme event type."""
    if extreme_type in NEGATIVE_TYPES:
        return -1.0, -1.0
    if extreme_type in POSITIVE_TYPES:
        return 1.0, 1.0
    raise ValueError("extreme_type 仅支持 NEGATIVE_TYPES 或 POSITIVE_TYPES")


def identify_extremes(
    climate_daily,
    *,
    extreme_type: str,
    scale: int = 30,
    nsp: int | bool = 30,
    si_method: SiMethod = "nonparametric",
    group_by: Literal["sum", "mean"] = "mean",
    pp_method: Literal["gringorten", "weibull"] = "gringorten",
    start_th: float | None = None,
    end_th: float | None = None,
    events_number_control: float = 0.1,
    events_days_control: float = 120.0,
    p: float = 0.05,
    distributions: Sequence[str] | None = None,
    coerce_365: bool = True,
) -> dict[str, Any]:
    """
    Single climate variable extreme event identification (SI → thresholds → PRM events → event summary).

    Parameters:
    - climate_daily: Daily-scale input (Series/DataFrame/xarray), preferably with DatetimeIndex.
    - extreme_type: Event type, must belong to NEGATIVE_TYPES/POSITIVE_TYPES in PRM_extreme_identification.
    - group_by: "sum" for precipitation, "mean" for temperature.
    - start_th / end_th: Pre-identification thresholds, auto-set by type if None (negative: -1/-1, positive: 1/1).

    Returns:
    {
        "climate_365", "si",
        "thresholds",
        "extreme_daily",
        "extreme_events",
    }
    """
    if si_method not in ("nonparametric", "best_distribution"):
        raise ValueError("si_method only supports 'nonparametric' or 'best_distribution'")

    # Auto thresholds
    if start_th is None or end_th is None:
        auto_start, auto_end = _default_threshold(extreme_type)
        start_th = auto_start if start_th is None else start_th
        end_th = auto_end if end_th is None else end_th

    # Normalize to 365/366 days
    climate_proc = to_365(climate_daily, how=group_by) if coerce_365 else climate_daily

    # Compute SI
    if si_method == "nonparametric":
        si = SI_nonparametric(
            climate_proc,
            scale=scale,
            group_by=group_by,
            NSP=nsp,
            pp_method=pp_method,
        )
    else:
        si = SI_best_distribution(
            climate_proc,
            scale=scale,
            group_by=group_by,
            NSP=nsp,
            distributions=distributions,
        )

    # Threshold search
    thresholds = remo_merg(
        extreme_type,
        si,
        scale=scale,
        start_th=start_th,
        end_th=end_th,
        events_number_control=events_number_control,
        events_days_control=events_days_control,
        p=p,
        distributions=distributions,
    )

    removal_raw, merging_raw = _get_thresholds(thresholds)
    removal, merging = _sanitize_thresholds(removal_raw, merging_raw, default_merging=-np.inf)

    # PRM identification + event table
    _, extreme_daily = PRM_extreme_identification(
        extreme_type,
        si,
        start_th,
        end_th,
        removal,
        merging,
    )
    extreme_events = daily_2_events(extreme_daily, extreme_type)

    return {
        "climate_365": climate_proc,
        "si": si,
        "thresholds": thresholds,
        "extreme_daily": extreme_daily,
        "extreme_events": extreme_events,
    }


def identify_compound_events(
    pre_daily,
    tem_daily,
    *,
    scale_pre: int = 30,
    scale_tem: int = 3,
    nsp: int | bool = 30,
    si_method: SiMethod = "nonparametric",
    group_by: Literal["sum", "mean"] = "mean",
    pp_method: Literal["gringorten", "weibull"] = "gringorten",
    start_th_d: float = -1.0,
    end_th_d: float = -1.0,
    start_th_h: float = 1.0,
    end_th_h: float = 1.0,
    events_number_control: float = 0.1,
    events_days_control: float = 120.0,
    p: float = 0.05,
    distributions: Sequence[str] | None = None,
    compound_type: int | str = 1,
    coerce_365: bool = True,
) -> dict[str, Any]:
    """
    Complete workflow from raw daily data to compound drought-heatwave event identification.

    Core Parameters:
    - pre_daily: Daily precipitation (Series/DataFrame/xarray), preferably with DatetimeIndex.
    - tem_daily: Daily mean temperature (Series/DataFrame/xarray).
    - si_method: "nonparametric" (default) or "best_distribution".
    - scale_pre / scale_tem: Accumulation window for SPI/SHI, default 30/3 days.
    - nsp: Non-stationary window (years), False means stationary.
    - start/end_th_d/h: Pre-identification thresholds for drought/heatwave.
    - compound_type: 1|2|3|4 or string alias, passed to identify_compound.

    Returns:
    {
        "pre_365", "tem_365",
        "si_pre", "si_tem",
        "drought_thresholds", "heatwave_thresholds",
        "drought_daily", "heatwave_daily",
        "drought_events", "heatwave_events",
        "compound_flags", "compound_daily",
    }
    """
    if si_method not in ("nonparametric", "best_distribution"):
        raise ValueError("si_method only supports 'nonparametric' or 'best_distribution'")

    pre_proc = to_365(pre_daily, how="sum") if coerce_365 else pre_daily
    tem_proc = to_365(tem_daily, how="mean") if coerce_365 else tem_daily

    if si_method == "nonparametric":
        si_pre = SI_nonparametric(
            pre_proc,
            scale=scale_pre,
            group_by="sum",
            NSP=nsp,
            pp_method=pp_method,
        )
        si_tem = SI_nonparametric(
            tem_proc,
            scale=scale_tem,
            group_by="mean",
            NSP=nsp,
            pp_method=pp_method,
        )
    else:
        si_pre = SI_best_distribution(
            pre_proc,
            scale=scale_pre,
            group_by="sum",
            NSP=nsp,
            distributions=distributions,
        )
        si_tem = SI_best_distribution(
            tem_proc,
            scale=scale_tem,
            group_by="mean",
            NSP=nsp,
            distributions=distributions,
        )

    drought_thresholds = remo_merg(
        "d",
        si_pre,
        scale=scale_pre,
        start_th=start_th_d,
        end_th=end_th_d,
        events_number_control=events_number_control,
        events_days_control=events_days_control,
        p=p,
        distributions=distributions,
    )
    heatwave_thresholds = remo_merg(
        "h",
        si_tem,
        scale=scale_tem,
        start_th=start_th_h,
        end_th=end_th_h,
        events_number_control=events_number_control,
        events_days_control=events_days_control,
        p=p,
        distributions=distributions,
    )

    drought_removal_raw, drought_merging_raw = _get_thresholds(drought_thresholds)
    heatwave_removal_raw, heatwave_merging_raw = _get_thresholds(heatwave_thresholds)

    drought_removal, drought_merging = _sanitize_thresholds(drought_removal_raw, drought_merging_raw, default_merging=-np.inf)
    heatwave_removal, heatwave_merging = _sanitize_thresholds(heatwave_removal_raw, heatwave_merging_raw, default_merging=-np.inf)

    _, drought_full = PRM_extreme_identification(
        "d",
        si_pre,
        start_th_d,
        end_th_d,
        drought_removal,
        drought_merging,
    )
    _, heatwave_full = PRM_extreme_identification(
        "h",
        si_tem,
        start_th_h,
        end_th_h,
        heatwave_removal,
        heatwave_merging,
    )

    drought_events = daily_2_events(drought_full, "d")
    heatwave_events = daily_2_events(heatwave_full, "h")
    compound_flags, compound_daily = identify_compound(drought_full, heatwave_full, compound_type)

    return {
        "pre_365": pre_proc,
        "tem_365": tem_proc,
        "si_pre": si_pre,
        "si_tem": si_tem,
        "drought_thresholds": drought_thresholds,
        "heatwave_thresholds": heatwave_thresholds,
        "drought_daily": drought_full,
        "heatwave_daily": heatwave_full,
        "drought_events": drought_events,
        "heatwave_events": heatwave_events,
        "compound_flags": compound_flags,
        "compound_daily": compound_daily,
    }


__all__ = ["identify_compound_events", "identify_extremes"]

