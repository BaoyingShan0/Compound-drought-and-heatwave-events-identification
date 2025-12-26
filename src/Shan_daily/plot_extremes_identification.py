import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import grid


def plot_extremes_identification(
    extreme_daily_all: pd.DataFrame, # a given period, including SI, flag_merged, order_merged
    start_th: float,
    extreme_type: str| None = None,
    fig_size: tuple[float, float] | None = (10, 4),
    save_path: str | None = None,
):
    """
    Plot extreme event identification results.
    Parameters:
    - extreme_daily_all: a given period, including SI, flag_merged, order_merged
    - start_th: the start threshold
    - extreme_type: the type of the extreme events
    - fig_size: figure size
    - save_path: path to save the figure
    """
    
    extremes_type_dict = {
        "d": "Drought",
        "dr": "Drought",
        "drought": "Drought",
        "c": "Coldwave",
        "cw": "Coldwave",
        "coldwave": "Coldwave",
        "p": "Pluvial",
        "pluvial": "Pluvial",
        "wetspell": "Pluvial",
        "h": "Heatwave",
        "heatwave": "Heatwave",
        "hotspell": "Heatwave"
    }

    if extreme_type in {"d", "dr", "c", "cw", "drought", "coldwave"}:
        extreme_color = [1, 0.8, 0]
    elif extreme_type in {"p", "h", "pluvial", "heatwave", "wetspell", "hotspell"}:
        extreme_color = [0.6, 0.8, 0.5]
    else:
        extreme_color = "gray"
    # Bar plot days that are in extreme events
    extreme_daily_2 = extreme_daily_all.copy()

    plt.figure(figsize=fig_size)
    mask = extreme_daily_2["flag_merged"].fillna(0).astype(bool)
    extreme_daily_2.loc[~mask, "SI"] = np.nan
    ba = plt.bar(extreme_daily_2.index, extreme_daily_2["SI"], color=extreme_color)
    ba.BarWidth = 1
   
    plt.plot(extreme_daily_all.index, extreme_daily_all["SI"], color=[0.85,0.33,0.10], linewidth=0.8)
    plt.plot(extreme_daily_all.index, np.full(len(extreme_daily_all.index), start_th), color="r", linestyle="--", linewidth=0.5)
    plt.xlim(extreme_daily_all.index[0], extreme_daily_all.index[-1])
    plt.ylabel("Standardized Index")
    plt.xlabel("Date")
    grid(True)
    if extreme_type is not None:
        plt.title(extremes_type_dict[extreme_type] + " identification results")
    else:
        plt.title("Extreme identification results")
    plt.grid(True)

    if save_path is not None:
        plt.savefig(save_path)
    plt.show()
