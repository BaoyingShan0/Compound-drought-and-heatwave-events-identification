# Nonparametric method to calculate standardized indices
# From: Farahmand, Alireza, and Amir AghaKouchak. "A Generalized Framework for Deriving Nonparametric Standardized Drought Indicators." Advances in Water Resources 76 (February 1, 2015): 140–45.
# https://doi.org/10.1016/j.advwatres.2014.11.012.
#
# 2025.12.22, by Baoying Shan

# Input Data:
# climate_data: daily climate data, pandas dataframe or xarray
# scale: days, like 15, 30 days; accumulation period
# NSP: variable for non-stationarity. If NSP=30, consider the past 30 years as "normal" condition
# If NSP=False (no input), data is assumed stationary, using the entire time series as normal condition.

# Output:
# SI: daily standardized index value, pandas dataframe or xarray


"""
Nonparametric standardized index based on empirical distribution (Gringorten plotting position).
- Accumulate daily data over a given timescale (rolling sum), compute probability p=(i-0.44)/(n+0.12) from empirical distribution, then apply standard normal inverse function to get SI.
- NSP=False means stationary assumption: for each calendar day, use the full-sample empirical distribution; NSP=30 (etc) means for trending days, use a moving window of ~NSP years.
Compare within same calendar day.
"""

from __future__ import annotations

from typing import Iterable
from bisect import bisect_left, bisect_right, insort

import numpy as np
import pandas as pd
from pandas.api.types import is_datetime64_any_dtype
from scipy.stats import norm

try:
    import xarray as xr
except ImportError:  # xarray 非必选
    xr = None

try:  # 可选的更快滚动实现
    import bottleneck as bn
except ImportError:
    bn = None

try:  # 可选 numba 加速滑窗
    import numba as nb
except ImportError:
    nb = None

def _rolling_sum(series: pd.Series, scale: int, group_by='mean') -> pd.Series:
    """Accumulate over given timescale, full window required. group_by='mean' for daily average, 'sum' for daily accumulation."""
    if bn is not None:
        arr = pd.Series(series).to_numpy(dtype=float)
        res = bn.move_sum(arr, window=scale, min_count=scale)
        if group_by == 'mean':
            res = res / scale
        return pd.Series(res, index=series.index)
    # Without bottleneck, use pandas rolling
    s = pd.Series(series)
    res = s.rolling(scale, min_periods=scale).sum()
    if group_by == 'mean':
        res = res / scale
    return res


def _prob_from_rank(rank: pd.Series | float, n: int, method: str) -> pd.Series | float:
    """Compute probability based on plotting position formula."""
    method = method.lower()
    if method.startswith("gring"):
        return (rank - 0.44) / (n + 0.12)
    if method.startswith("weibull"):
        return rank / (n + 1)
    raise ValueError("pp_method only supports 'gringorten' or 'weibull'")


def _seasonal_global_si(acc: pd.Series, pp_method: str) -> pd.Series:
    """Stationary case: compute empirical distribution SI separately for each calendar day-of-year."""
    idx = pd.DatetimeIndex(acc.index)
    doy = _dayofyear_no_leap(idx)
    out = pd.Series(np.nan, index=acc.index)
    for day in np.unique(doy):
        mask = doy == day
        vals = acc[mask]
        nonnan = vals[vals.notna()].sort_values()
        n = len(nonnan)
        if n == 0:
            continue
        ranks = nonnan.rank(method="max")
        p = _prob_from_rank(ranks, n, pp_method)
        si = norm.ppf(p)
        out.loc[nonnan.index] = si
    return out


def _pp_method_code(pp_method: str) -> int:
    """Return internal code for plotting position formula."""
    method = pp_method.lower()
    if method.startswith("gring"):
        return 0
    if method.startswith("weibull"):
        return 1
    raise ValueError("pp_method only supports 'gringorten' or 'weibull'")


def _prob_from_rank_code(rank: int, n: int, code: int) -> float:
    if code == 0:  # gringorten
        return (rank - 0.44) / (n + 0.12)
    # weibull
    return rank / (n + 1.0)


def _local_probs_py(values: np.ndarray, window: int, code: int) -> np.ndarray:
    """Pure Python sliding window (sorted list maintenance) for probability, fallback for numba."""
    n = len(values)
    probs = np.full(n, np.nan, dtype=float)
    if window <= 0 or n < window:
        return probs

    sorted_vals: list[float] = []
    # Initialize first window
    for i in range(window):
        v = values[i]
        if np.isfinite(v):
            insort(sorted_vals, float(v))

    for i in range(n):
        if i < window:
            new = values[i]
            if np.isfinite(new):
                insort(sorted_vals, float(new))
        if i >= window:
            old = values[i - window]
            if np.isfinite(old):
                pos = bisect_left(sorted_vals, float(old))
                if pos < len(sorted_vals) and abs(sorted_vals[pos] - old) < 1e-12:
                    sorted_vals.pop(pos)
            new = values[i]
            if np.isfinite(new):
                insort(sorted_vals, float(new))

        val = values[i]
        if not (np.isfinite(val)):
            continue
        m = len(sorted_vals)
        if m == 0:
            continue
        rank = bisect_right(sorted_vals, float(val))
        probs[i] = _prob_from_rank_code(rank, m, code)
    return probs


if nb is not None:
    @nb.njit
    def _insert_sorted(lst, val):
        i = 0
        n = len(lst)
        while i < n and lst[i] <= val:
            i += 1
        lst.insert(i, val)

    @nb.njit
    def _remove_first(lst, val):
        n = len(lst)
        for i in range(n):
            if abs(lst[i] - val) < 1e-12:
                lst.pop(i)
                return

    @nb.njit
    def _searchsorted_right(lst, val):
        lo = 0
        hi = len(lst)
        while lo < hi:
            mid = (lo + hi) // 2
            if val < lst[mid]:
                hi = mid
            else:
                lo = mid + 1
        return lo

    @nb.njit
    def _local_probs_nb(values: np.ndarray, window: int, code: int) -> np.ndarray:
        n = len(values)
        probs = np.empty(n, dtype=np.float64)
        probs[:] = np.nan
        if window <= 0 or n < window:
            return probs

        lst = nb.typed.List.empty_list(nb.types.float64)
        for i in range(window):
            v = values[i]
            if np.isfinite(v):
                _insert_sorted(lst, v)
        
        for i in range(n):
            if i<window:
                new = values[i]
                if np.isfinite(new):
                    _insert_sorted(lst, new)
                    
            if i >= window:
                old = values[i - window]
                if np.isfinite(old):
                    _remove_first(lst, old)
                new = values[i]
                if np.isfinite(new):
                    _insert_sorted(lst, new)

            val = values[i]
            if not (np.isfinite(val)):
                continue
            m = len(lst)
            if m == 0:
                continue
            rank = _searchsorted_right(lst, val)
            if code == 0:
                probs[i] = (rank - 0.44) / (m + 0.12)  # gringorten
            else:
                probs[i] = rank / (m + 1.0)  # weibull
        

        return probs


def _seasonal_local_si(
    acc: pd.Series, window_years: int, pp_method: str, target_days: set[int] | None = None
) -> pd.Series:
    """For specified days, compute SI using recent years' sliding window (by year samples); return NaN for others."""
    idx = pd.DatetimeIndex(acc.index)
    doy = _dayofyear_no_leap(idx)
    values = pd.Series(acc).to_numpy(dtype=float)
    code = _pp_method_code(pp_method)

    probs = np.full(len(values), np.nan, dtype=float)
    unique_days = np.unique(doy) if target_days is None else np.array(sorted(target_days))
    for day in unique_days:
        if target_days is not None and day not in target_days:
            continue
        pos = np.nonzero(doy == day)[0]
        if len(pos) == 0:
            continue
        day_vals = values[pos]
        if nb is not None:
            day_probs = _local_probs_nb(day_vals, window_years, code)
        else:
            day_probs = _local_probs_py(day_vals, window_years, code)
        for k, p in zip(pos, day_probs):
            probs[k] = p

    si_vals = norm.ppf(probs)
    return pd.Series(si_vals, index=acc.index)


def _to_series(data: Iterable) -> pd.Series:
    if isinstance(data, pd.Series):
        return data
    return pd.Series(data)


def _trend_days_by_mk(acc: pd.Series, alpha: float = 0.05, min_samples: int = 8) -> set[int]:
    """Perform MK test for each day-of-year, return set of day indices with significant trend."""
    try:
        import pymannkendall as mk  # type: ignore[import-not-found]
    except ImportError as exc:  # only required when needed
        raise ImportError("Need to install pymannkendall for non-stationary testing") from exc

    idx = pd.DatetimeIndex(acc.index)
    doy = _dayofyear_no_leap(idx)
    trend_days: set[int] = set()

    for day in np.unique(doy):
        vals = acc[(doy == day) & np.isfinite(acc)]
        if len(vals) < min_samples:
            continue
        res = mk.original_test(vals.values)
        if getattr(res, "p", 1.0) <= alpha and getattr(res, "trend", "no trend") != "no trend":
            trend_days.add(day)
    return trend_days

def _dayofyear_no_leap(dt_index: pd.DatetimeIndex) -> pd.Series:
    """Compute day-of-year excluding leap year offset."""
    doy = dt_index.dayofyear
    return doy - ((dt_index.is_leap_year) & (dt_index.month > 2))

def SI_nonparametric(
    climate_data, 
    scale: int, 
    group_by: Literal["mean", "sum"] = "mean", 
    NSP: bool | int = False, 
    pp_method: Literal["gringorten", "weibull"] = "gringorten"
    ) -> pd.Series | pd.DataFrame | xr.DataArray | xr.Dataset:
    """
    Compute nonparametric standardized index.
    Parameters:
    - climate_data: Daily-scale data (pandas Series/single or multi-column DataFrame/iterable convertible to Series, or xarray DataArray/Dataset with time dimension).
    - scale: Accumulation window (days), using rolling sum or rolling mean.
    - group_by: Grouping method, 'mean' for daily average, 'sum' for daily accumulation.
    - NSP: Whether to perform non-stationary distribution calculation. False means stationary (use full-sample empirical distribution for each calendar day);
      if integer (years), first perform MK trend test for each calendar day (alpha=0.05), use sliding window of ~NSP years for trending days, full-sample otherwise.
    - pp_method: Plotting position formula; "gringorten" (default) or "weibull" (rank/(n+1)).
    Returns: pandas Series/DataFrame or xarray with same index as input.
    """
    # pandas DataFrame: compute by column and return DataFrame (skip date columns)
    if isinstance(climate_data, pd.DataFrame):
        cols = [
            col
            for col in climate_data.columns
            if str(col).lower() != "date" and not is_datetime64_any_dtype(climate_data[col])
        ]
        if not cols:  # only date column, return empty frame
            return pd.DataFrame(index=climate_data.index)
        if len(cols) > 1:
            out = {col: SI_nonparametric(climate_data[col], scale, group_by, NSP, pp_method) for col in cols}
            return pd.DataFrame(out, index=climate_data.index)
        # single column falls through to generic Series branch
        climate_data = climate_data[cols[0]]

    # xarray DataArray
    if xr is not None and isinstance(climate_data, xr.DataArray):
        if "time" not in climate_data.dims:
            raise ValueError("DataArray must contain time dimension")
        series = pd.Series(climate_data.values, index=pd.DatetimeIndex(climate_data["time"].values))
        res = SI_nonparametric(series, scale, group_by, NSP, pp_method)
        return xr.DataArray(res.values, coords={"time": climate_data["time"]}, dims=["time"], name=climate_data.name)

    # xarray Dataset: process each variable (only 1D time supported)
    if xr is not None and isinstance(climate_data, xr.Dataset):
        out = {}
        for name, da in climate_data.data_vars.items():
            if ("time" not in da.dims) or (len(da.dims) != 1):
                raise ValueError("Dataset only supports 1D variables with time dimension only")
            out[name] = SI_nonparametric(da, scale, group_by, NSP, pp_method)
        return xr.Dataset(out, coords={"time": climate_data["time"]})

    series = _to_series(climate_data)
    acc = _rolling_sum(series, scale, group_by)

    if NSP is False:
        return _seasonal_global_si(acc, pp_method)

    if not isinstance(NSP, (int, float)) or NSP <= 0:
        raise ValueError("NSP should be a positive number of years or False")

    # Need datetime index for day-wise MK test
    try:
        doy = _dayofyear_no_leap(pd.DatetimeIndex(acc.index))
    except Exception as exc:
        raise ValueError("Datetime index required for non-stationary detection") from exc

    window_years = int(round(NSP))

    trend_days = _trend_days_by_mk(acc, alpha=0.05)
    # If no significant trend, directly return full-sample SI (by day)
    if not trend_days:
        return _seasonal_global_si(acc, pp_method)

    global_si = _seasonal_global_si(acc, pp_method)
    local_si = _seasonal_local_si(acc, window_years, pp_method, target_days=trend_days)

    trend_mask = doy.isin(trend_days)
    
    out = global_si.copy()
    out.loc[trend_mask] = local_si.loc[trend_mask]
    return out
