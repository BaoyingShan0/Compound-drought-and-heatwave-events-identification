from __future__ import annotations

import numpy as np
import pandas as pd


NEGATIVE_TYPES = {"d", "c", "dr", "cw", "drought", "coldwave"}


def daily_2_events(extreme_daily: pd.DataFrame, extreme_type: str) -> pd.DataFrame:
    """
    Summarize daily-scale identification results (PRM_extreme_identification's full DataFrame) into event table.

    Output columns:
    - num: event number (starting from 1)
    - duration: number of days
    - severity: SI sum (inverted for drought/cold types)
    - intensity: SI average (inverted for drought/cold types)
    - neighbor_duration: length of non-event period between current event and previous event
    - inter_arrival: duration + neighbor_duration
    """
    if not isinstance(extreme_daily, pd.DataFrame):
        raise ValueError("extreme_daily should be a DataFrame")
    if "order_merged" not in extreme_daily or "flag_merged" not in extreme_daily or "SI" not in extreme_daily:
        raise ValueError("extreme_daily must contain SI/flag_merged/order_merged columns")

    # Only keep event rows (flag_merged==1 and order is negative)
    order = extreme_daily["order_merged"].to_numpy()
    flag = extreme_daily["flag_merged"].to_numpy()
    si = extreme_daily["SI"].to_numpy(dtype=float)

    event_mask = (flag == 1) & np.isfinite(order) & (order < 0)
    if not event_mask.any():
        return pd.DataFrame(columns=["num", "duration", "severity", "intensity", "neighbor_duration", "inter_arrival"])

    unique_orders = np.array(sorted(np.unique(order[event_mask]), reverse=True), dtype=float)  # negative values descending: -1 ... -n
    event_nums = {o: i + 1 for i, o in enumerate(unique_orders)}

    rows = []
    for o in unique_orders:
        mask = order == o
        idx_positions = np.flatnonzero(mask)
        duration = len(idx_positions)
        severity = np.nansum(si[mask])
        intensity = severity / duration if duration > 0 else np.nan
        rows.append(
            {
                "num": event_nums[o],
                "duration": duration,
                "severity": severity,
                "intensity": intensity,
                "_first_idx": idx_positions[0],
                "_last_idx": idx_positions[-1],
            }
        )

    events_df = pd.DataFrame(rows).set_index("num")

    # Compute neighbor intervals
    neighbor_durations = [np.nan]
    inter_arrivals = [np.nan]
    for i in range(2, len(events_df) + 1):
        cur_first = events_df.loc[i, "_first_idx"]
        prev_last = events_df.loc[i - 1, "_last_idx"]
        gap = cur_first - prev_last - 1
        neighbor_durations.append(gap)
        inter_arrivals.append(gap + events_df.loc[i, "duration"])

    events_df["neighbor_duration"] = neighbor_durations
    events_df["inter_arrival"] = inter_arrivals

    # Invert (drought/cold)
    if extreme_type in NEGATIVE_TYPES:
        events_df["severity"] = -events_df["severity"]
        events_df["intensity"] = -events_df["intensity"]

    events_df = events_df[["duration", "severity", "intensity", "neighbor_duration", "inter_arrival"]]
    events_df.insert(0, "num", events_df.index)
    events_df.reset_index(drop=True, inplace=True)
    return events_df