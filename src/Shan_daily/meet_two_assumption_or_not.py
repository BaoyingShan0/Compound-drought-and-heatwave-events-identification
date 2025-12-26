from __future__ import annotations
from typing import Sequence
import numpy as np
import pandas as pd
from scipy import stats

"""
检验极端事件序列是否满足两个假设：
1) 到达间隔服从指数分布（泊松过程）；
2) 事件严重度（severity）服从 GEV 分布。
若两者都通过，再挑选最优分布（AIC 最小）。

输入:
- extreme_daily: DataFrame，需包含 SI、flag_merged、order_merged 列；order_merged<0 表示事件编号。
- extreme_type: "d"/"c" 等负向事件需要将 severity 取反，其余视为正向。
- p: KS 检验显著性水平，默认 0.05。
- distributions: 可选，候选分布名称列表（对应 scipy.stats），为空时使用默认集合。

输出:
- assumption1(bool): 间隔是否服从指数分布
- assumption2(bool): 严重度是否服从 GEV
- best_dist(str|None): 若两者均通过，返回最优分布名称，否则 None
"""

NEGATIVE_TYPES = {"d", "dr", "c", "cw", "drought", "coldwave"}


def _fit_best_distribution(
    values: np.ndarray,
    dist_candidates: Sequence[tuple[str, stats.rv_continuous]],
) -> tuple[str, stats.rv_continuous, tuple] | None:
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return None

    best: tuple[str, stats.rv_continuous, tuple] | None = None
    best_aic = np.inf
    for name, dist in dist_candidates:
        try:
            params = dist.fit(finite)
            k = len(params)
            logpdf = dist.logpdf(finite, *params)
            if not np.all(np.isfinite(logpdf)):
                continue
            ll = float(np.sum(logpdf))
            aic = 2 * k - 2 * ll
        except Exception:
            continue
        if aic < best_aic:
            best_aic = aic
            best = (name, dist, params)
    return best



# step 1: first calculate the severity of extreme spells and inter-arrival time, based
# on the order_merged and flag_merged columns
def _build_events(extreme_daily: pd.DataFrame, extreme_type: str) -> pd.DataFrame:
    if not isinstance(extreme_daily, pd.DataFrame):
        raise ValueError("extreme_daily 应为 DataFrame")
    required = {"order_merged", "flag_merged", "SI"}
    if not required.issubset(extreme_daily.columns):
        raise ValueError("extreme_daily 需包含 SI/flag_merged/order_merged 列")

    order = extreme_daily["order_merged"].to_numpy()
    flag = extreme_daily["flag_merged"].to_numpy()
    si = extreme_daily["SI"].to_numpy(dtype=float)

    event_mask = (flag == 1) & np.isfinite(order) & (order < 0)
    if not event_mask.any():
        return pd.DataFrame(columns=["severity", "inter_arrival"])

    unique_orders = np.array(sorted(np.unique(order[event_mask]), reverse=True), dtype=float)  # -1, -2, ...
    rows = []
    for num, o in enumerate(unique_orders, start=1):
        mask = order == o
        idx_positions = np.flatnonzero(mask)
        duration = len(idx_positions)
        severity = np.nansum(si[mask])
        rows.append(
            {
                "num": num,
                "duration": duration,
                "severity": severity,
                "_first_idx": idx_positions[0],
                "_last_idx": idx_positions[-1],
            }
        )

    events_df = pd.DataFrame(rows).set_index("num")

    inter_arrivals = [np.nan]
    for i in range(2, len(events_df) + 1):
        cur_first = events_df.loc[i, "_first_idx"]
        prev_last = events_df.loc[i - 1, "_last_idx"]
        gap = cur_first - prev_last - 1
        inter_arrivals.append(gap + events_df.loc[i, "duration"])
    events_df["inter_arrival"] = inter_arrivals

    if extreme_type in NEGATIVE_TYPES:
        events_df["severity"] = -events_df["severity"]

    return events_df[["severity", "inter_arrival"]]


def _ks_test(series: pd.Series, dist: stats.rv_continuous, alpha: float) -> tuple[bool, float, tuple]:
    finite = series[np.isfinite(series)]
    if len(finite) == 0:
        return False, np.nan, ()
    params = dist.fit(finite)
    cdf = lambda x: dist.cdf(x, *params)
    stat = stats.kstest(finite, cdf)
    return stat.pvalue >= alpha, stat.pvalue, params


def meet_two_assumption_or_not(
    extreme_daily: pd.DataFrame,
    extreme_type: str,
    p: float = 0.05,
    distributions: Sequence[tuple[str, stats.rv_continuous]] | None = None,
):
    """
    返回 (assumptions_passed, best_dist)：
    - assumptions_passed: 间隔服从指数分布且严重度服从 GEV
    - best_dist: 若两者都通过且提供 distributions，则返回 AIC 最优分布名称，否则 None
    """
    events_df = _build_events(extreme_daily, extreme_type)
    if events_df.empty:
        return False, None

    assumption1, _, _ = _ks_test(events_df["inter_arrival"], stats.expon, p)

    assumption2 = False
    best_dist_name: str | None = None
    if assumption1:
        assumption2, _, _ = _ks_test(events_df["severity"], stats.genextreme, p)
        if assumption2 and distributions is not None:
            # distributions 已经是 list[tuple[str, rv_continuous]] 格式
            best = _fit_best_distribution(events_df["severity"].to_numpy(), distributions)
            if best is not None:
                best_dist_name = best[0]

    return assumption1 and assumption2, best_dist_name