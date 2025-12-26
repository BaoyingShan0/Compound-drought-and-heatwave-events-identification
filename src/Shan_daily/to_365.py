"""
Normalize daily climate data to 365 days:
- Input supports pandas Series/DataFrame (requires datetime index or date column) or xarray DataArray/Dataset (requires time dimension).
- how="sum" adds Feb 29 to Feb 28 (precipitation and other accumulations); how="mean" averages the two days (temperature and other averages).
- After processing each leap year, remove Feb 29. Output type matches input, pandas preserves original date column or DatetimeIndex.
"""

from __future__ import annotations

from typing import Literal, Tuple

import pandas as pd

try:
    import xarray as xr
except ImportError:  # xarray 非必选
    xr = None


How = Literal["sum", "mean"]


def _to_datetime_index(df: pd.DataFrame | pd.Series, date_col: str | None) -> Tuple[pd.DataFrame | pd.Series, bool]:
    """Ensure datetime index, return (object, whether to restore date column)."""
    if isinstance(df, pd.Series):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Series requires DatetimeIndex")
        return df, False

    needs_restore = False
    if isinstance(df.index, pd.DatetimeIndex):
        return df, False

    col = date_col or "date"
    if col not in df.columns:
        raise ValueError(f"Date column {col} not found")

    df = df.copy()
    df.index = pd.to_datetime(df[col])
    df.drop(columns=[col], inplace=True)
    needs_restore = True
    return df, needs_restore


def _combine_feb29_pd(obj: pd.DataFrame | pd.Series, how: How, date_col: str | None) -> pd.DataFrame | pd.Series:
    obj, restore_col = _to_datetime_index(obj, date_col)
    op = (lambda a, b: (a + b) / 2) if how == "mean" else (lambda a, b: a + b)

    idx = obj.index
    feb29 = (idx.month == 2) & (idx.day == 29)
    if not feb29.any():
        return _restore(obj, restore_col)

    years = pd.Index(idx[feb29].year).unique()
    obj = obj.copy()

    for year in years:
        m29 = feb29 & (idx.year == year)
        m28 = (idx.month == 2) & (idx.day == 28) & (idx.year == year)
        if not m29.any():
            continue
        if not m28.any():
            # If Feb 28 missing, directly delete Feb 29
            obj = obj[~m29]
            idx = obj.index
            feb29 = (idx.month == 2) & (idx.day == 29)
            continue

        if isinstance(obj, pd.Series):
            v28 = obj.loc[m28].iloc[0]
            v29 = obj.loc[m29].iloc[0]
            obj.loc[m28] = op(v28, v29)
        else:
            num_cols = obj.select_dtypes(include="number").columns
            v28 = obj.loc[m28, num_cols]
            v29 = obj.loc[m29, num_cols]
            obj.loc[m28, num_cols] = op(v28.to_numpy(), v29.to_numpy())

        obj = obj[~m29]
        idx = obj.index
        feb29 = (idx.month == 2) & (idx.day == 29)

    return _restore(obj, restore_col)


def _restore(obj: pd.DataFrame | pd.Series, restore_col: bool) -> pd.DataFrame | pd.Series:
    """Restore date column (if not originally index)."""
    if restore_col and isinstance(obj, pd.DataFrame):
        obj = obj.copy()
        obj.insert(0, "date", obj.index)
    return obj


def to_365(data, how: How = "sum", date_col: str | None = None):
    """
    Compress leap year daily data to 365 days.

    Parameters
    - data: pandas Series/DataFrame (DatetimeIndex or date column), or xarray DataArray/Dataset (with time).
    - how: "sum" | "mean", determines how to merge Feb 28 and Feb 29.
    - date_col: When DataFrame doesn't have DatetimeIndex, specify date column name (default "date").
    """
    if how not in ("sum", "mean"):
        raise ValueError('how must be "sum" or "mean"')

    # pandas
    if isinstance(data, (pd.Series, pd.DataFrame)):
        return _combine_feb29_pd(data, how, date_col)

    # xarray
    if xr is not None and isinstance(data, xr.DataArray):
        processed = _combine_feb29_pd(data.to_pandas(), how, date_col=None)
        return xr.DataArray(processed, dims=["time"], name=data.name)

    if xr is not None and isinstance(data, xr.Dataset):
        processed = _combine_feb29_pd(data.to_dataframe(), how, date_col=None)
        return xr.Dataset.from_dataframe(processed)

    raise TypeError("Only supports pandas Series/DataFrame or xarray DataArray/Dataset")