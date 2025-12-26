from __future__ import annotations

from typing import Sequence
import numpy as np
import pandas as pd
from scipy import stats
from .PRM_extreme_identification import PRM_extreme_identification
from .daily_2_events import daily_2_events
from .proximity import proximity
from .meet_two_assumption_or_not import meet_two_assumption_or_not

"""
## Important note: 
# 2025.12.24, by Baoying Shan
# Compared to the removal and merging algorithm in HESS paper, 
# we only use min AIC to select the best distribution for severity among the distributions

Algorithm:
Generate candidate space:
Run PRM_extreme_identification with initial thresholds once → get event duration distribution.
removal candidates: duration+1, truncate to ≤ 5*scale, add 0 (no removal).
For each removal candidate, compute inter-event "proximity" → get merging candidates, filter overly large values (>3*removal or >5*scale; if removal=0, only >5*scale), add -inf (no merging).
Form matrix REMO_d (rows) and MERG_d (row-column combinations).

Grid search + pruning:
Traverse undecided cells in matrix (removal, merging combinations).
For each combination, run PRM_extreme_identification once to get event sequence.
Compute annual event count and annual event days; if control conditions not met, prune entire row/column to 0 (infeasible).
If control conditions met, use meet_two_assumption_or_not to check "inter-arrival~exponential and severity~GEV" assumptions; if passed, mark cell as 0.5 and prune lower-right submatrix (monotonicity pruning for larger event counts); if failed, mark as 0.
Select optimal:
Collect all combinations satisfying both assumptions, select by maximum "annual event count" (return its removal, merging, event count, event days). Return NaN if no feasible solution.
"""
#%%

# give the searching space for the removal and merging thresholds
def _get_searching_space(
    extreme_type: str,
    SI_values: pd.Series,
    scale: int,
    start_th: float,
    end_th: float,
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """
    Generate search space for removal and merging thresholds.
    Returns:
    - REMO_d: possible removal threshold array (including 0 for no removal)
    - MERG_D: dict, key is row index (aligned with REMO_d order), value is optional merging threshold array for that removal (including -inf for no merging)
    """
    _, extreme_daily = PRM_extreme_identification(extreme_type, SI_values, start_th, end_th, -np.inf, -np.inf)
    events = daily_2_events(extreme_daily, extreme_type)
    duration = events["duration"]  # duration is the number of days in the spell

    remo_candidates = np.unique(np.sort(duration + 1))
    remo_candidates = remo_candidates[remo_candidates <= scale * 5]
    remo_candidates = np.concatenate([[0], remo_candidates])  # add 0, meaning no removal

    pairs: list[tuple[float, np.ndarray]] = []
    for remo_d in remo_candidates:
        prox = proximity(extreme_type, SI_values, start_th, end_th, remo_d, -np.inf)
        if prox is None or len(prox) == 0:
            continue
        merge_d = np.unique(np.sort(np.ceil(prox)))
        if remo_d != 0:
            cc = (merge_d > 3 * remo_d) | (merge_d > 5 * scale)
        else:
            cc = merge_d > 5 * scale
        merge_d = merge_d[~cc]
        if len(merge_d) == 0:
            continue
        merge_d = np.concatenate([[-np.inf], merge_d])  # add -inf, meaning no merging
        pairs.append((remo_d, merge_d))

    if not pairs:
        return np.array([]), {}

    REMO_d = np.array([p[0] for p in pairs])
    MERG_D = {i: p[1] for i, p in enumerate(pairs)}
    return REMO_d, MERG_D

def _default_candidates() -> list[tuple[str, stats.rv_continuous]]:
    return [
        ("norm", stats.norm),
        ("expon", stats.expon),
        ("gamma", stats.gamma),
        ("genextreme", stats.genextreme), # GEV
        ("invgauss", stats.invgauss),
        ("logistic", stats.logistic),
        ("fisk", stats.fisk),  # log-logistic
        ("lognorm", stats.lognorm),
        ("burr", stats.burr),
        ("gumbel_r", stats.gumbel_r),
    ]
#%%
def remo_merg(
    extreme_type: str,
    SI_values: pd.Series,
    scale: int,
    start_th: float,
    end_th: float,
    events_number_control: float = 0.1,
    events_days_control: float = 120,
    p: float = 0.05,
    distributions: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Search for optimal removal and merging threshold combination, return DataFrame with optimal combination.
    """
    dist_candidates = [(dist, getattr(stats, dist)) for dist in distributions] if distributions is not None else _default_candidates()

    REMO_d, MERG_D = _get_searching_space(extreme_type, SI_values, scale, start_th, end_th)
    if len(REMO_d) == 0 or len(MERG_D) == 0:
        return pd.DataFrame(
            {
                "removal_threshold": [np.nan],
                "merging_threshold": [np.nan],
                "events_per_year": [np.nan],
                "days_per_year": [np.nan],
            }
        )

    MERG_d = np.full((len(REMO_d), max(len(v) for v in MERG_D.values())), np.nan)
    for idx, merge_d in MERG_D.items():
        MERG_d[idx, : len(merge_d)] = merge_d

    ID_R = len(REMO_d)
    ID_M = MERG_d.shape[1]
    mat_id = np.full((ID_R, ID_M), np.nan)  # 0 means impossible, 0.5 means possible, 1 means best
    mat_id[np.isnan(MERG_d)] = 0

    times = 0
    id_r = 0  # starting search point for removal threshold
    id_m = 0
    years = int(SI_values.index.year.nunique()) if hasattr(SI_values.index, "year") else len(SI_values) / 365
    Best_dist_dr = np.full((ID_R, ID_M), None, dtype=object)  # store best distribution name (string)
    Assumption_all = np.full((ID_R, ID_M), np.nan)  # whether two assumptions hold

    while np.any(np.isnan(mat_id)):
        if id_r >= ID_R or id_m >= ID_M:
            break

        times += 1
        remo_d = REMO_d[id_r]
        merg_d = MERG_d[id_r, id_m]

        _, extreme_daily_all = PRM_extreme_identification(extreme_type, SI_values, start_th, end_th, remo_d, merg_d)
        events_per_yr = -float(np.min(extreme_daily_all["order_merged"])) / years
        days_per_yr = float(np.sum(extreme_daily_all["flag_merged"])) / years

        if events_per_yr < events_number_control:
            mat_id[id_r:, id_m] = 0
        elif days_per_yr > events_days_control:
            mat_id[id_r, id_m:] = 0
        else:
            assumption, best_dist_dr = meet_two_assumption_or_not(
                extreme_daily_all, 
                extreme_type, 
                p, 
                distributions=dist_candidates)
            Best_dist_dr[id_r, id_m] = best_dist_dr
            Assumption_all[id_r, id_m] = assumption

            if assumption and best_dist_dr == "genextreme":
                mat_id[id_r:, id_m:] = 0
                mat_id[id_r, id_m] = 1
                print(f"remo_d: {remo_d}, merg_d: {merg_d}, best_dist_dr: {best_dist_dr}, assumption: {assumption}")
            elif assumption:
                mat_id[id_r, id_m] = 0.5
                #print(f"remo_d: {remo_d}, merg_d: {merg_d}, best_dist_dr: {best_dist_dr}, assumption: {assumption}")
            else:
                mat_id[id_r, id_m] = 0

        # choose next position
        while id_r < ID_R and id_m < ID_M and not np.isnan(mat_id[id_r, id_m]):
            col_nan = np.isnan(mat_id[:, id_m])
            if col_nan.any():
                # prefer smallest id_r with nan in current column
                nan_rows = np.where(col_nan)[0]
                if (nan_rows < id_r).any():
                    id_r = nan_rows[nan_rows < id_r][-1]
                elif (nan_rows > id_r).any():
                    id_r = nan_rows[nan_rows > id_r][0]
                else:
                    id_r = id_r  # no change
            else:
                id_m += 1
                id_r = 0
            if id_m >= ID_M:
                break

    id_r_best, id_m_best = np.where((Assumption_all == 1) & (Best_dist_dr == "genextreme"))
    if id_r_best.size == 0:
        id_r_best, id_m_best = np.where(Assumption_all == 1)
    if id_r_best.size == 0:
        return pd.DataFrame(
            {
                "removal_threshold": [np.nan],
                "merging_threshold": [np.nan],
                "events_per_year": [np.nan],
                "days_per_year": [np.nan],
            }
        )

    events_per_yr = np.full(id_r_best.size, np.nan)
    days_per_yr = np.full(id_r_best.size, np.nan)
    for i in range(id_r_best.size):
        _, extreme_daily_all = PRM_extreme_identification(
            extreme_type,
            SI_values,
            start_th,
            end_th,
            REMO_d[id_r_best[i]],
            MERG_d[id_r_best[i], id_m_best[i]],
        )
        events_per_yr[i] = -float(np.min(extreme_daily_all["order_merged"])) / years
        days_per_yr[i] = float(np.sum(extreme_daily_all["flag_merged"])) / years

    ind = int(np.argmax(events_per_yr))
    best_combination = pd.DataFrame(
        {
            "removal_threshold": [REMO_d[id_r_best[ind]]],
            "merging_threshold": [MERG_d[id_r_best[ind], id_m_best[ind]]],
            "events_per_year": [events_per_yr[ind]],
            "days_per_year": [days_per_yr[ind]],
        }
    )
    return best_combination


