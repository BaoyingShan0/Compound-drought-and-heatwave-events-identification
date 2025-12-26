from __future__ import annotations

import numpy as np
import pandas as pd

try:
    import xarray as xr
except ImportError:  # 可选依赖
    xr = None

from .PRM_extreme_identification import PRM_extreme_identification

NEGATIVE_TYPES = {"d", "c", "dr", "cw", "drought", "coldwave"}
POSITIVE_TYPES = {"p", "h", "pluvial", "heatwave", "wetspell", "hotspell"}
def _event_ids_from_orders(orders: np.ndarray, flags: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """将 order_merged/flag_merged 转为正整数事件编号及对应的 order 值。"""
    mask = flags == 1
    if not mask.any():
        return np.array([], dtype=int), np.array([], dtype=float)
    unique_orders = np.unique(orders[mask])
    # order=0,-1,-2,... -> 1,2,3,...（保持顺序）
    order_to_id = {o: i + 1 for i, o in enumerate(sorted(unique_orders))}
    ids = np.array([order_to_id[o] for o in unique_orders], dtype=int)
    return ids, unique_orders


def _proximity_one_series(series: pd.Series, extreme_type: str, start_th: float, end_th: float, REMO, MERG):
    """单列 proximity 计算，返回事件接近度 Series（index=事件编号，从1开始）。"""
    if not isinstance(series.index, pd.DatetimeIndex):
        raise ValueError("series 需为带 DatetimeIndex 的 pandas Series")

    # 调用 PRM 获取事件掩码与序号（PRM 内部处理符号）
    flag_merged, full = PRM_extreme_identification(
        extreme_type, series, start_th, end_th, REMO=REMO, MERG=MERG
    )

    orders = full["order_merged"].to_numpy(dtype=float)
    flags_not_in_spells = (1-full["flag_merged"]).to_numpy(dtype=float)

    # 与 PRM 内部一致的阈值与指数符号，用于 proximity 求和
    index_eff = series.to_numpy(dtype=float)
    start_th_eff = start_th
    if extreme_type in NEGATIVE_TYPES:
        index_eff = -index_eff
        start_th_eff = -start_th
    elif extreme_type in POSITIVE_TYPES:
        pass
    else:
        raise ValueError("extreme_type 不支持")

    neighbor_ids, neighbor_orders = _event_ids_from_orders(orders, flags_not_in_spells)
    if len(neighbor_ids) == 0:
        return pd.Series(dtype=float)

    # proximity = sum(start_th_eff - index_eff) over each事件
    prox = np.full(neighbor_ids.max(), np.nan, dtype=float)
    for eid, o in zip(neighbor_ids, neighbor_orders):
        mask = orders == o
        prox[eid - 1] = np.nansum(start_th_eff - index_eff[mask])

    return pd.Series(prox, index=pd.Index(range(1, len(prox) + 1), name="neighbor_id"))


def proximity(extreme_type, SI_values, start_th, end_th, REMO=0, MERG=-np.inf):
    """
    MATLAB proximity.m 的 Python 版。

    输入支持 pandas Series / DataFrame（忽略日期列）/ xarray DataArray 或 Dataset（需含 time）。
    返回：
    - 输入为 Series 或 DataArray：返回 proximity Series（index 为事件编号，从 1 开始）。
    - 输入为多列 DataFrame：返回 proximity DataFrame，列与原列名一致。
    - 输入为 xarray Dataset：返回 Dataset，各变量为原变量名 + "_proximity"。
    """
    # xarray DataArray
    if xr is not None and isinstance(SI_values, xr.DataArray):
        if "time" not in SI_values.dims:
            raise ValueError("DataArray 需包含 time 维度")
        series = pd.Series(SI_values.values, index=pd.DatetimeIndex(SI_values["time"].values))
        return _proximity_one_series(series, extreme_type, start_th, end_th, REMO, MERG)

    # xarray Dataset
    if xr is not None and isinstance(SI_values, xr.Dataset):
        out = {}
        neighbor_ids = []
        for name, da in SI_values.data_vars.items():
            if ("time" not in da.dims) or (len(da.dims) != 1):
                raise ValueError("Dataset 仅支持一维 time 变量")
            series = pd.Series(da.values, index=pd.DatetimeIndex(da["time"].values))
            prox = _proximity_one_series(series, extreme_type, start_th, end_th, REMO, MERG)
            if len(prox) > 0:
                out[f"{name}_proximity"] = ("neighbor_id", prox.values)
                neighbor_ids.append(prox.index)
        # 使用最长的 neighbor_id 作为坐标
        if neighbor_ids:
            max_len = max(len(idx) for idx in neighbor_ids)
            neighbor_coord = pd.Index(range(1, max_len + 1), name="neighbor_id")
        else:
            neighbor_coord = pd.Index([], name="neighbor_id", dtype=int)
        return xr.Dataset(out, coords={"neighbor_id": neighbor_coord})

    # pandas DataFrame
    if isinstance(SI_values, pd.DataFrame):
        cols = [
            col
            for col in SI_values.columns
            if str(col).lower() != "date" and not pd.api.types.is_datetime64_any_dtype(SI_values[col])
        ]
        if not cols:
            return pd.DataFrame(index=pd.Index([], name="event_id"))
        if len(cols) > 1:
            prox_cols = {}
            max_len = 0
            for col in cols:
                prox_series = _proximity_one_series(SI_values[col], extreme_type, start_th, end_th, REMO, MERG)
                prox_cols[col] = prox_series
                max_len = max(max_len, len(prox_series))
            # 对齐长度
            for col, s in prox_cols.items():
                if len(s) < max_len:
                    prox_cols[col] = s.reindex(range(1, max_len + 1))
            return pd.DataFrame(prox_cols, index=pd.Index(range(1, max_len + 1), name="event_id"))
        # 单列
        return _proximity_one_series(SI_values[cols[0]], extreme_type, start_th, end_th, REMO, MERG)

    # pandas Series
    if isinstance(SI_values, pd.Series):
        return _proximity_one_series(SI_values, extreme_type, start_th, end_th, REMO, MERG)

    # 其他可迭代 -> Series（需外部保证 index 为 datetime）
    series = pd.Series(SI_values)
    if not isinstance(series.index, pd.DatetimeIndex):
        raise ValueError("SI_values 需提供 DatetimeIndex")
    return _proximity_one_series(series, extreme_type, start_th, end_th, REMO, MERG)

