from __future__ import annotations

import numpy as np
import pandas as pd

TYPE_MAP = {
    1: "ex1-and-ex2",
    2: "ex1-or-ex2",
    3: "ex1-cond-ex2",
    4: "ex2-cond-ex1",
}
TYPE_ALIAS = {
    "and": "ex1-and-ex2",
    "intersection": "ex1-and-ex2",
    "union": "ex1-or-ex2",
    "or": "ex1-or-ex2",
    "ex1": "ex1-cond-ex2",
    "ex2": "ex2-cond-ex1",
}
TYPE_NAMES = set(TYPE_MAP.values())

__all__ = ["identify_compound"]


def _normalize_type(type_value: int | str) -> str:
    """Map numeric/alias to one of four combination types."""
    if isinstance(type_value, str):
        name = TYPE_ALIAS.get(type_value.strip().lower(), type_value.strip().lower())
        if name not in TYPE_NAMES:
            raise ValueError("type only supports 1/2/3/4 or corresponding strings.")
        return name
    if isinstance(type_value, (int, np.integer)):
        if type_value not in TYPE_MAP:
            raise ValueError("type only supports 1/2/3/4 or corresponding strings.")
        return TYPE_MAP[int(type_value)]
    raise ValueError("type only supports 1/2/3/4 or corresponding strings.")


def _prepare_series(df: pd.DataFrame, name: str) -> tuple[np.ndarray, np.ndarray]:
    """Extract SI and flag_merged from PRM_extreme_identification output."""
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame.")

    columns = {str(c).lower(): c for c in df.columns}
    si_col = columns.get("si")
    flag_col = columns.get("flag_merged") or columns.get("flag")
    if si_col is None or flag_col is None:
        raise ValueError(f"{name} must contain columns 'SI' and 'flag_merged'.")

    si = df[si_col].to_numpy(dtype=float)
    flags = df[flag_col].to_numpy(dtype=float)
    if len(si) != len(flags):
        raise ValueError(f"{name}'s SI and flag_merged must have the same length.")
    return si, flags


def _label_runs(flags: np.ndarray) -> np.ndarray:
    """
    Event segments: 1,2,3...; non-event segments: 0,-1,-2... (incrementing/decrementing order).
    """
    orders = np.zeros_like(flags, dtype=float)
    aa = 0.0  # non-event segment ID (decreasing)
    bb = 0.0  # event segment ID (increasing)

    for i, f in enumerate(flags):
        if i == 0:
            orders[i] = bb if f == 1 else aa
            continue
        prev = flags[i - 1]
        if f == 1 and prev == 1:
            orders[i] = bb
        elif f == 1 and prev == 0:
            bb += 1
            orders[i] = bb
        elif f == 0 and prev == 1:
            aa -= 1
            orders[i] = aa
        else:
            orders[i] = aa
    return orders


def _build_event_flags(Index1_0: np.ndarray, Index2_0: np.ndarray, type_name: str) -> np.ndarray:
    """Generate 0/1 event flags based on four combination types."""
    if type_name == "ex1-and-ex2":
        return (
            (~np.isnan(Index1_0))
            & (~np.isnan(Index2_0))
            & (Index1_0 != 0)
            & (Index2_0 != 0)
        ).astype(float)

    if type_name == "ex1-or-ex2":
        return (
            ~((np.isnan(Index1_0) & np.isnan(Index2_0)) | ((Index1_0 == 0) & (Index2_0 == 0)))
        ).astype(float)

    if type_name == "ex1-cond-ex2":
        return (~np.isnan(Index1_0) & (Index1_0 != 0)).astype(float)

    if type_name == "ex2-cond-ex1":
        return (~np.isnan(Index2_0) & (Index2_0 != 0)).astype(float)

    raise ValueError("Unknown combination type, use 1/2/3/4 or corresponding strings.")


def identify_compound(
    ex1_daily: pd.DataFrame,
    ex2_daily: pd.DataFrame,
    type: int | str = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Input: Output from PRM_extreme_identification (containing SI, flag_merged).
    Returns: (compound_df, compound_daily).
    """
    if len(ex1_daily) != len(ex2_daily):
        raise ValueError("ex1_daily and ex2_daily must have the same length.")

    type_name = _normalize_type(type)
    si1, flag1 = _prepare_series(ex1_daily, "ex1_daily")
    si2, flag2 = _prepare_series(ex2_daily, "ex2_daily")

    Index1_0 = si1.copy()
    Index2_0 = si2.copy()
    Index1_0[flag1 == 0] = np.nan
    Index2_0[flag2 == 0] = np.nan

    flags = _build_event_flags(Index1_0, Index2_0, type_name)
    orders = _label_runs(flags)

    compound = np.zeros((len(flags), 2), dtype=float)
    for order_id in np.unique(orders[orders > 0]):
        mask = orders == order_id
        diff = Index1_0[mask] - Index2_0[mask]
        has_full = np.all(~np.isnan(Index1_0))
        mean_close = not np.isnan(diff).all() and (np.abs(np.nanmean(diff)) < 1)
        has_overlap = not np.isnan(diff).all()
        if (has_full and mean_close) or (not has_full and has_overlap):
            compound[mask, 0] = 1
            compound[mask, 1] = compound[:, 1].max() + 1

    compound_df = pd.DataFrame(
        {
            "is_compound": compound[:, 0].astype(int),
            "compound_order": compound[:, 1].astype(int),
        },
        index=ex1_daily.index,
    )

    mask_compound = compound_df["is_compound"].to_numpy(dtype=bool)
    si1_compound = np.where(mask_compound, si1, 0.0)
    si2_compound = np.where(mask_compound, si2, 0.0)

    compound_daily = pd.DataFrame(
        {
            "date": ex1_daily.index if isinstance(ex1_daily.index, pd.DatetimeIndex) else ex1_daily.index,
            "ex1_SI": si1_compound,
            "ex2_SI": si2_compound,
            "is_compound": compound_df["is_compound"].to_numpy(),
            "compound_order": compound_df["compound_order"].to_numpy(),
        },
        index=ex1_daily.index,
    )

    return compound_df, compound_daily

