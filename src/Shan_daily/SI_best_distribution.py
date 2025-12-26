from __future__ import annotations

"""
Compute standardized indices based on best-fit distribution (selected by AIC).
Similar to SI_nonparametric, but fits 10 common distributions to each day's data,
selects the one with minimum AIC, then uses its CDF through standard normal inverse to get SI.

Candidate distributions: normal, exponential, gamma, GEV(genextreme),
inverse Gaussian(invgauss), logistic, log-logistic(fisk), log-normal(lognorm),
Burr, extreme value(gumbel_r).
"""

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from pandas.api.types import is_datetime64_any_dtype
from scipy import stats
from scipy.stats import norm

try:
    import xarray as xr
except ImportError:  # 可选依赖
    xr = None

try:  # 可选更快 rolling
    import bottleneck as bn  # type: ignore[import-not-found]
except ImportError:
    bn = None


def _rolling_sum(series: pd.Series, scale: int, group_by: str = "mean") -> pd.Series:
    """Accumulate over given timescale, full window required. group_by='mean' for daily average, 'sum' for daily accumulation."""
    s = pd.Series(series)
    if bn is not None:
        arr = s.to_numpy(dtype=float)
        res = bn.move_sum(arr, window=scale, min_count=scale)
        if group_by == "mean":
            res = res / scale
        return pd.Series(res, index=s.index)
    rolled = s.rolling(scale, min_periods=scale).sum()
    if group_by == "mean":
        rolled = rolled / scale
    return rolled


def _to_series(data: Iterable) -> pd.Series:
    if isinstance(data, pd.Series):
        return data
    return pd.Series(data)


def _dayofyear_no_leap(dt_index: pd.DatetimeIndex) -> pd.Series:
    """Compute day-of-year excluding leap year offset."""
    doy = dt_index.dayofyear
    return doy - ((dt_index.is_leap_year) & (dt_index.month > 2))


def _fit_best_distribution(
    values: np.ndarray,
    dist_candidates: Sequence[tuple[str, stats.rv_continuous]],
) -> tuple[str, stats.rv_continuous, tuple] | None:
    """Fit candidate distributions to given sample, select best by AIC, return (name, dist, params)."""
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


def _cdf_to_si(values: np.ndarray, dist: stats.rv_continuous, params: tuple) -> np.ndarray:
    """Convert to standard normal SI using distribution CDF."""
    p = dist.cdf(values, *params)
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return norm.ppf(p)


def _seasonal_global_si(
    acc: pd.Series,
    dist_candidates: Sequence[tuple[str, stats.rv_continuous]],
    min_samples: int,
) -> pd.Series:
    """Stationary case: fit best distribution and compute SI for each calendar day."""
    idx = pd.DatetimeIndex(acc.index)
    doy = _dayofyear_no_leap(idx)
    out = pd.Series(np.nan, index=acc.index, dtype=float)

    for day in np.unique(doy):
        mask = doy == day
        vals = acc[mask]
        finite = vals[vals.notna()]
        if len(finite) < min_samples:
            continue

        all_nonneg = np.all(finite.to_numpy(dtype=float) >= 0)
        zeros_cnt = int((finite == 0).sum())
        if all_nonneg and zeros_cnt > 0:
            q = zeros_cnt / len(finite)
            nonzero = finite[finite > 0]

            # First handle zero values' probability
            zero_prob = float(np.clip(q, 1e-12, 1 - 1e-12))
            if zeros_cnt > 0:
                out.loc[finite[finite == 0].index] = norm.ppf(zero_prob)

            # For non-zero values, use conditional probability q + (1-q)*F
            if len(nonzero) >= min_samples:
                best = _fit_best_distribution(nonzero.to_numpy(dtype=float), dist_candidates)
            else:
                best = None
            if best is None:
                continue
            _, dist, params = best
            nz_vals = vals[(vals.notna()) & (vals.to_numpy(dtype=float) > 0)]
            if len(nz_vals) > 0:
                p = dist.cdf(nz_vals.to_numpy(dtype=float), *params)
                p = q + (1 - q) * p
                p = np.clip(p, 1e-12, 1 - 1e-12)
                out.loc[nz_vals.index] = norm.ppf(p)
            continue

        best = _fit_best_distribution(finite.to_numpy(dtype=float), dist_candidates)
        if best is None:
            continue
        _, dist, params = best
        si_vals = _cdf_to_si(vals.to_numpy(dtype=float), dist, params)
        out.loc[vals.index] = si_vals
    return out


def _seasonal_local_si(
    acc: pd.Series,
    window_years: int,
    dist_candidates: Sequence[tuple[str, stats.rv_continuous]],
    min_samples: int,
    target_days: set[int] | None = None,
) -> pd.Series:
    """For specified days, fit distribution using ~window_years samples and compute SI; return NaN for others."""
    idx = pd.DatetimeIndex(acc.index)
    doy = _dayofyear_no_leap(idx)
    out = pd.Series(np.nan, index=acc.index, dtype=float)

    unique_days = np.unique(doy) if target_days is None else np.array(sorted(target_days))
    for day in unique_days:
        if target_days is not None and day not in target_days:
            continue
        pos = np.nonzero(doy == day)[0]
        if len(pos) == 0:
            continue
        day_vals = acc.iloc[pos].to_numpy(dtype=float)

        for j, k in enumerate(pos):
            start = max(0, j - window_years + 1)
            window_vals = day_vals[start : j + 1]
            finite = window_vals[np.isfinite(window_vals)]
            if len(finite) < min_samples:
                continue

            all_nonneg = np.all(finite >= 0)
            zeros_cnt = int(np.sum(finite == 0))
            val = day_vals[j]
            # When all valid values are non-negative, assign probability q=z/M for zero values (z=zeros, M=sample size), 
            # use conditional probability q+(1−q)F for non-zero values (F is best-fit distribution CDF for non-zero samples), 
            # then apply standard normal inverse to get SI.
            if all_nonneg and zeros_cnt > 0:
                q = zeros_cnt / len(finite)
                nonzero = finite[finite > 0]

                if np.isfinite(val) and abs(val) < 1e-12:
                    p = q
                else:
                    if len(nonzero) >= min_samples:
                        best = _fit_best_distribution(nonzero, dist_candidates)
                    else:
                        best = None
                    if best is None:
                        continue
                    _, dist, params = best
                    F = dist.cdf(val, *params)
                    if not np.isfinite(F):
                        continue
                    p = q + (1 - q) * F
            else:
                best = _fit_best_distribution(finite, dist_candidates)
                if best is None:
                    continue
                _, dist, params = best
                p = dist.cdf(val, *params)
                if not np.isfinite(p):
                    continue

            p = float(np.clip(p, 1e-12, 1 - 1e-12))
            out.iloc[k] = norm.ppf(p)
    return out


def _trend_days_by_mk(acc: pd.Series, alpha: float = 0.05, min_samples: int = 8) -> set[int]:
    """Perform MK test for each day-of-year, return set of day indices with significant trend."""
    try:
        import pymannkendall as mk  # type: ignore[import-not-found]
    except ImportError as exc:
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


def _default_candidates() -> list[tuple[str, stats.rv_continuous]]:
    """Return default candidate distribution list."""
    return [
        ("norm", stats.norm),
        ("expon", stats.expon),
        ("gamma", stats.gamma),
        ("genextreme", stats.genextreme),
        ("invgauss", stats.invgauss),
        ("logistic", stats.logistic),
        ("fisk", stats.fisk),  # log-logistic
        ("lognorm", stats.lognorm),
        ("burr", stats.burr),
        ("gumbel_r", stats.gumbel_r),  # EV
    ]


def SI_best_distribution(
    climate_data,
    scale: int,
    group_by: str = "mean",
    NSP: bool | int = False,
    distributions: Sequence[str] | None = None,
    min_samples: int = 5,
):
    """
    
    Compute standardized index based on best-fit distribution (selected by AIC).

    Parameters:
    - climate_data: Daily-scale data (pandas Series/single or multi-column DataFrame/iterable convertible to Series,
      or xarray DataArray/Dataset with time dimension).
    - scale: Accumulation window (days), using rolling sum or rolling mean.
    - group_by: 'mean' (default, daily average) or 'sum' (daily accumulation).
    - NSP: Whether to perform non-stationary distribution calculation. False means stationary (fit best distribution using full sample for each calendar day);
      if positive integer, use sliding window of ~NSP years for trending days, full sample otherwise.
    - distributions: Candidate distribution list (name, scipy.stats distribution object), defaults to 10 common distributions.
    - min_samples: Minimum sample size for fitting distributions, default 5.

    Returns: pandas Series/DataFrame or xarray with same index as input.
    """
    dist_candidates = [(dist, getattr(stats, dist)) for dist in distributions] if distributions is not None else _default_candidates()

    # pandas DataFrame: compute by column and return DataFrame (skip date columns)
    if isinstance(climate_data, pd.DataFrame):
        cols = [
            col
            for col in climate_data.columns
            if str(col).lower() != "date" and not is_datetime64_any_dtype(climate_data[col])
        ]
        if not cols:
            return pd.DataFrame(index=climate_data.index)
        if len(cols) > 1:
            out = {
                col: SI_best_distribution(climate_data[col], scale, group_by, NSP, dist_candidates, min_samples)
                for col in cols
            }
            return pd.DataFrame(out, index=climate_data.index)
        climate_data = climate_data[cols[0]]

    # xarray DataArray
    if xr is not None and isinstance(climate_data, xr.DataArray):
        if "time" not in climate_data.dims:
            raise ValueError("DataArray must contain time dimension")
        series = pd.Series(climate_data.values, index=pd.DatetimeIndex(climate_data["time"].values))
        res = SI_best_distribution(series, scale, group_by, NSP, dist_candidates, min_samples)
        return xr.DataArray(res.values, coords={"time": climate_data["time"]}, dims=["time"], name=climate_data.name)

    # xarray Dataset: process each variable (only 1D time supported)
    if xr is not None and isinstance(climate_data, xr.Dataset):
        out = {}
        for name, da in climate_data.data_vars.items():
            if ("time" not in da.dims) or (len(da.dims) != 1):
                raise ValueError("Dataset only supports 1D variables with time dimension only")
            out[name] = SI_best_distribution(da, scale, group_by, NSP, dist_candidates, min_samples)
        return xr.Dataset(out, coords={"time": climate_data["time"]})

    series = _to_series(climate_data)
    acc = _rolling_sum(series, scale, group_by)

    if NSP is False:
        return _seasonal_global_si(acc, dist_candidates, min_samples)

    if not isinstance(NSP, (int, float)) or NSP <= 0:
        raise ValueError("NSP should be a positive number of years or False")

    try:
        doy = _dayofyear_no_leap(pd.DatetimeIndex(acc.index))
    except Exception as exc:
        raise ValueError("Datetime index required for non-stationary detection") from exc

    window_years = int(round(NSP))
    trend_days = _trend_days_by_mk(acc, alpha=0.05)
    if not trend_days:
        return _seasonal_global_si(acc, dist_candidates, min_samples)

    global_si = _seasonal_global_si(acc, dist_candidates, min_samples)
    local_si = _seasonal_local_si(acc, window_years, dist_candidates, min_samples, target_days=trend_days)

    trend_mask = doy.isin(trend_days)
    out = global_si.copy()
    out.loc[trend_mask] = local_si.loc[trend_mask]
    return out



