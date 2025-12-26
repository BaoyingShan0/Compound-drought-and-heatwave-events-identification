from __future__ import annotations

import numpy as np
import pandas as pd

try:
    import xarray as xr
except ImportError:  # 可选依赖
    xr = None

NEGATIVE_TYPES = {"d", "c", "dr", "cw", "drought", "coldwave"}
POSITIVE_TYPES = {"p", "h", "pluvial", "heatwave", "wetspell", "hotspell"}
def _compute_orders(flags: np.ndarray) -> np.ndarray:
    """Compute phase numbers from flag sequence (dry periods as non-positive numbers, wet periods as non-negative)."""
    orders = np.empty_like(flags, dtype=float)
    aa = 0.0  # non-dry period ID (increasing)
    bb = 0.0  # dry period ID (decreasing)
    for i, f in enumerate(flags):
        if i == 0:
            orders[i] = bb if f == 1 else aa
            continue
        prev = flags[i - 1]
        if f == 1 and prev == 1:
            orders[i] = bb
        elif f == 1 and prev == 0:
            bb -= 1
            orders[i] = bb
        elif f == 0 and prev == 1:
            aa += 1
            orders[i] = aa
        else:
            orders[i] = aa
    return orders


def _process_one_series(
    extreme_type: str,
    series: pd.Series,
    start_th: float,
    end_th: float,
    REMO: float,
    MERG: float,
) -> pd.DataFrame:
    """Perform PRM extreme event identification on a single Series."""
    original_si = series.to_numpy(dtype=float)
    index = original_si.copy()

    if extreme_type in NEGATIVE_TYPES:
        index = -index
        start_th = -start_th
        end_th = -end_th
    elif extreme_type in POSITIVE_TYPES:
        pass
    else:
        raise ValueError("extreme_type only supports NEGATIVE_TYPES or POSITIVE_TYPES")

    n = len(index)

    # STEP 1: Pre-identification
    flag_pre = np.zeros(n, dtype=float)
    for i, val in enumerate(index):
        if i == 0:
            flag_pre[i] = 1.0 if val >= start_th else 0.0
        else:
            if val >= start_th:
                flag_pre[i] = 1.0
            elif flag_pre[i - 1] == 1.0 and val > end_th:
                flag_pre[i] = 1.0
            else:
                flag_pre[i] = 0.0
    order_pre = _compute_orders(flag_pre)

    # STEP 2: Remove short dry periods (if REMO provided)
    flag_remove = np.zeros(n, dtype=float)
    if REMO is not None:
        max_order = int(np.nanmax(order_pre)) if n > 0 else 0
        if np.isnan(max_order):
            max_order = 0
        for i in range(1, max_order + 1):
            mask = order_pre == -i
            if not mask.any():
                continue
            flag_remove[mask] = 0.0 if mask.sum() < REMO else 1.0
    order_remove = _compute_orders(flag_remove)

    # STEP 3: Merge adjacent events (if MERG provided)
    flag_merge = flag_remove.copy()
    if MERG is not None:
        max_order_rm = int(np.nanmax(order_remove)) if n > 0 else 0
        if np.isnan(max_order_rm):
            max_order_rm = 0
        for i in range(1, max_order_rm + 1):
            mask = order_remove == i
            if not mask.any():
                continue
            if np.nansum(start_th - index[mask]) < MERG:
                flag_merge[mask] = 1.0
            else:
                flag_merge[mask] = 0.0
    order_merge = _compute_orders(flag_merge)

    # STEP 4: Restore original SI, sync NaN
    nan_mask = np.isnan(original_si)
    flag_pre[nan_mask] = np.nan
    flag_remove[nan_mask] = np.nan
    flag_merge[nan_mask] = np.nan

    return pd.DataFrame(
        {
            "SI": original_si,
            "flag_pre": flag_pre,
            "order_pre": order_pre,
            "flag_removed": flag_remove,
            "order_removed": order_remove,
            "flag_merged": flag_merge,
            "order_merged": order_merge,
        },
        index=series.index,
    )


def _extract_flag_merged(result):
    """Extract primary return value flag_merged from full result, preserving input type consistency."""
    if isinstance(result, pd.DataFrame):
        if isinstance(result.columns, pd.MultiIndex):
            flag_cols = result.xs("flag_merged", axis=1, level=1, drop_level=False)
            if flag_cols.shape[1] == 1:
                flag_cols = flag_cols.droplevel(0, axis=1)
            return flag_cols
        return result["flag_merged"]

    if xr is not None and isinstance(result, xr.Dataset):
        flag_vars = {k: v for k, v in result.data_vars.items() if k.endswith("flag_merged")}
        return xr.Dataset(flag_vars, coords=result.coords)

    return result


def PRM_extreme_identification(
    extreme_type: str,
    SI_values,
    start_th: float,
    end_th: float,
    REMO: float | None = 0,
    MERG: float | None = -np.inf,
):
    """
    Replicate MATLAB PRM_extreme_identification, supporting Series / DataFrame /
    iterable / xarray DataArray or Dataset (must contain time dimension).
    Returns two objects: (flag_merged subset, full result).
    - Single-column DataFrame: flag_merged as Series.
    - Multi-column DataFrame: flag_merged as DataFrame with only flag_merged (preserving multi-level columns).
    - Input xarray DataArray: first return Dataset with only flag_merged, second return full Dataset.
    - Input xarray Dataset: first return Dataset with only *_flag_merged variables, second return full Dataset.
    """
    # xarray DataArray
    if xr is not None and isinstance(SI_values, xr.DataArray):
        if "time" not in SI_values.dims:
            raise ValueError("DataArray must contain time dimension")
        series = pd.Series(SI_values.values, index=pd.DatetimeIndex(SI_values["time"].values))
        df = _process_one_series(extreme_type, series, start_th, end_th, REMO, MERG)
        ds = xr.Dataset({col: (("time",), df[col].values) for col in df.columns}, coords={"time": SI_values["time"]})
        return _extract_flag_merged(ds), ds

    # xarray Dataset
    if xr is not None and isinstance(SI_values, xr.Dataset):
        out = {}
        time_coord = SI_values["time"]
        for name, da in SI_values.data_vars.items():
            if ("time" not in da.dims) or (len(da.dims) != 1):
                raise ValueError("Dataset only supports 1D time variables")
            series = pd.Series(da.values, index=pd.DatetimeIndex(da["time"].values))
            df = _process_one_series(extreme_type, series, start_th, end_th, REMO, MERG)
            for col in df.columns:
                out[f"{name}_{col}"] = (("time",), df[col].values)
        ds = xr.Dataset(out, coords={"time": time_coord})
        return _extract_flag_merged(ds), ds

    # pandas DataFrame
    if isinstance(SI_values, pd.DataFrame):
        cols = [
            col
            for col in SI_values.columns
            if str(col).lower() != "date" and not pd.api.types.is_datetime64_any_dtype(SI_values[col])
        ]
        if not cols:
            empty = pd.DataFrame(index=SI_values.index)
            return _extract_flag_merged(empty), empty
        if len(cols) > 1:
            pieces = {
                col: _process_one_series(extreme_type, SI_values[col], start_th, end_th, REMO, MERG)
                for col in cols
            }
            full = pd.concat(pieces, axis=1)
            return _extract_flag_merged(full), full
        # single column
        full = _process_one_series(extreme_type, SI_values[cols[0]], start_th, end_th, REMO, MERG)
        return _extract_flag_merged(full), full

    # pandas Series
    if isinstance(SI_values, pd.Series):
        if not isinstance(SI_values.index, pd.DatetimeIndex):
            raise ValueError("SI_values must be a pandas Series with DatetimeIndex")
        full = _process_one_series(extreme_type, SI_values, start_th, end_th, REMO, MERG)
        return _extract_flag_merged(full), full

    # other iterable -> Series (external caller must ensure index is datetime)
    series = pd.Series(SI_values)
    if not isinstance(series.index, pd.DatetimeIndex):
        raise ValueError("SI_values must provide DatetimeIndex")
    full = _process_one_series(extreme_type, series, start_th, end_th, REMO, MERG)
    return _extract_flag_merged(full), full