"""
Demonstration of all functions in Shan_daily package for extreme event identification.
"""

from pathlib import Path

import pandas as pd

# Main entry functions
from Shan_daily import identify_compound_events, identify_extremes

# Individual components for advanced usage
from Shan_daily import (
    to_365,
    SI_nonparametric,
    SI_best_distribution,
    PRM_extreme_identification,
    remo_merg,
    daily_2_events,
    identify_compound,
    plot_extremes_identification,
    identify_compound_events,
    identify_extremes,
)


def main():
    """Run comprehensive demo of Shan_daily package."""
    
    # ========================================
    # 1. Load example data
    # ========================================
    print("=" * 60)
    print("1. Loading example data...")
    print("=" * 60)
    
    base = Path(__file__).resolve().parent
    pre_path = base / "data/daily_pre.csv"
    tem_path = base / "data/tem_daily.csv"
    
    pre_raw = pd.read_csv(pre_path)
    tem_raw = pd.read_csv(tem_path)
    
    print(f"Precipitation data shape: {pre_raw.shape}")
    print(f"Temperature data shape: {tem_raw.shape}")
    print()
    
    # ========================================
    # 2. Main Entry Function: Compound Events
    # ========================================
    print("=" * 60)
    print("2. Running identify_compound_events (main workflow)...")
    print("=" * 60)
    
    result_compound = identify_compound_events(
        pre_daily=pre_raw,
        tem_daily=tem_raw,
        scale_pre=30,           # 30-day SPI
        scale_tem=3,            # 3-day SHI
        nsp=30,                 # 30-year non-stationary window
        si_method="nonparametric",  # or "best_distribution"
        start_th_d=-1.0,        # drought start threshold
        end_th_d=-1.0,          # drought end threshold
        start_th_h=1.0,         # heatwave start threshold
        end_th_h=1.0,           # heatwave end threshold
        compound_type=1,        # 1: AND, 2: OR, 3: drought-cond-heatwave, 4: heatwave-cond-drought
        coerce_365=True,        # normalize to 365 days
    )
    
    print("Compound events workflow completed!")
    print(f"Keys in result: {list(result_compound.keys())}")
    print(f"Number of drought events: {len(result_compound['drought_events'])}")
    print(f"Number of heatwave events: {len(result_compound['heatwave_events'])}")
    print(f"Number of compound days: {result_compound['compound_daily']['is_compound'].sum()}")
    print()
    
    # Display sample results
    print("Sample drought events:")
    print(result_compound['drought_events'].head())
    print()
    
    print("Sample compound days:")
    print(result_compound['compound_daily'][result_compound['compound_daily']['is_compound'] == 1].head())
    print()
    
    # ========================================
    # 3. Main Entry Function: Single Variable Extremes
    # ========================================
    print("=" * 60)
    print("3. Running identify_extremes (single variable)...")
    print("=" * 60)
    
    # Identify drought events only
    result_drought = identify_extremes(
        climate_daily=pre_raw,
        extreme_type="d",       # "d" for drought
        scale=30,
        nsp=30,
        si_method="nonparametric",
        group_by="sum",         # "sum" for precipitation
        start_th=-1.0,
        end_th=-1.0,
        coerce_365=True,
    )
    
    print("Single variable extreme identification completed!")
    print(f"Keys in result: {list(result_drought.keys())}")
    print(f"Number of drought events: {len(result_drought['extreme_events'])}")
    print()
    
    # ========================================
    # 4. Advanced: Individual Component Functions
    # ========================================
    print("=" * 60)
    print("4. Testing individual component functions...")
    print("=" * 60)
    
    # 4.1 to_365: Normalize to 365 days
    print("\n4.1 to_365: Normalizing leap year data...")
    pre_365 = to_365(pre_raw, how="sum")
    tem_365 = to_365(tem_raw, how="mean")
    print(f"Pre-normalization shape: {pre_raw.shape}, Post-normalization: {pre_365.shape}")
    
    # 4.2 SI_nonparametric: Compute standardized indices
    print("\n4.2 SI_nonparametric: Computing standardized indices...")
    spi = SI_nonparametric(pre_365['pre'], scale=30, group_by='sum', NSP=30)
    shi = SI_nonparametric(tem_365['tem_mean'], scale=3, group_by='mean', NSP=30)
    print(f"SPI range: [{spi.min():.2f}, {spi.max():.2f}]")
    print(f"SHI range: [{shi.min():.2f}, {shi.max():.2f}]")
    
    # 4.3 SI_best_distribution: Alternative SI method
    print("\n4.3 SI_best_distribution: Computing SI with best-fit distribution...")
    distributions = ["norm", "gamma", "genextreme", "burr"]
    spi_best = SI_best_distribution(
        pre_365['pre'], 
        scale=30, 
        group_by='sum', 
        NSP=30,
        distributions=distributions
    )
    print(f"SPI (best dist) range: [{spi_best.min():.2f}, {spi_best.max():.2f}]")
    
    # 4.4 remo_merg: Search for optimal thresholds
    print("\n4.4 remo_merg: Searching for optimal removal/merging thresholds...")
    drought_thresholds = remo_merg(
        "d",
        spi,
        scale=30,
        start_th=-1.0,
        end_th=-1.0,
        events_number_control=0.1,
        events_days_control=120.0,
        p=0.05,
    )
    print("Optimal thresholds found:")
    print(drought_thresholds)
    
    # 4.5 PRM_extreme_identification: Identify extreme events
    print("\n4.5 PRM_extreme_identification: Identifying extreme events with PRM...")
    removal_th = drought_thresholds.iloc[0]['removal_threshold']
    merging_th = drought_thresholds.iloc[0]['merging_threshold']
    
    _, drought_daily = PRM_extreme_identification(
        "d",
        spi,
        start_th=-1.0,
        end_th=-1.0,
        REMO=removal_th,
        MERG=merging_th,
    )
    print(f"Drought daily flags shape: {drought_daily.shape}")
    print(f"Number of drought days: {drought_daily['flag_merged'].sum()}")
    
    # 4.6 daily_2_events: Convert daily flags to event summary
    print("\n4.6 daily_2_events: Converting daily data to event summary...")
    drought_events = daily_2_events(drought_daily, "d")
    print(f"Number of drought events: {len(drought_events)}")
    print("Sample events:")
    print(drought_events.head())
    
    # 4.7 identify_compound: Identify compound events from two extremes
    print("\n4.7 identify_compound: Identifying compound events...")
    
    # Get heatwave daily flags
    _, heatwave_daily = PRM_extreme_identification(
        "h",
        shi,
        start_th=1.0,
        end_th=1.0,
        REMO=None,
        MERG=-float('inf'),
    )
    
    # Identify compound events (type 1: AND)
    compound_flags, compound_daily = identify_compound(
        drought_daily,
        heatwave_daily,
        type=1  # or "and", "intersection"
    )
    
    print(f"Number of compound days: {compound_flags['is_compound'].sum()}")
    print("Sample compound days:")
    print(compound_daily[compound_daily['is_compound'] == 1].head())
    
    # ========================================
    # 5. Visualization
    # ========================================
    print("\n" + "=" * 60)
    print("5. Creating visualizations...")
    print("=" * 60)
    
    #%% Plot first year of drought identification results
    plot_start = 1000
    plot_end = 1365
    
    print("\n5.1 Plotting drought identification results (first year)...")
    plot_extremes_identification(
        extreme_daily_all=drought_daily.iloc[plot_start:plot_end],
        start_th=-1.0,
        extreme_type="d",
        fig_size=(12, 3),
        #save_path=base / "outputs" / "drought_identification_year1.png",
    )
    
    print("\n5.2 Plotting heatwave identification results (first year)...")
    plot_extremes_identification(
        extreme_daily_all=heatwave_daily.iloc[plot_start:plot_end],
        start_th=1.0,
        extreme_type="h",
        fig_size=(12, 5),
        #save_path=base / "outputs" / "heatwave_identification_year1.png",
    )
    
    # ========================================
    # 6. Testing Different Compound Types
    # ========================================
    print("\n" + "=" * 60)
    print("6. Testing different compound event types...")
    print("=" * 60)
    
    compound_types = {
        1: "AND (intersection)",
        2: "OR (union)",
        3: "Drought-conditional-heatwave",
        4: "Heatwave-conditional-drought"
    }
    
    for ctype, description in compound_types.items():
        _, comp_daily = identify_compound(drought_daily, heatwave_daily, type=ctype)
        n_compound = comp_daily['is_compound'].sum()
        print(f"Type {ctype} ({description}): {n_compound} compound days")
    
    print()
    
    # ========================================
    # 7. Save Results
    # ========================================
    print("=" * 60)
    print("7. Saving results...")
    print("=" * 60)
    
    output_dir = base / "outputs"
    output_dir.mkdir(exist_ok=True)
    
    # Save event summaries
    result_compound['drought_events'].to_csv(output_dir / "drought_events.csv", index=False)
    result_compound['heatwave_events'].to_csv(output_dir / "heatwave_events.csv", index=False)
    result_compound['compound_daily'].to_csv(output_dir / "compound_daily.csv", index=False)
    
    # Save thresholds
    result_compound['drought_thresholds'].to_csv(output_dir / "drought_thresholds.csv", index=False)
    result_compound['heatwave_thresholds'].to_csv(output_dir / "heatwave_thresholds.csv", index=False)
    
    print("Results saved to:", output_dir)
    print()
    
    # ========================================
    # 8. Summary Statistics
    # ========================================
    print("=" * 60)
    print("8. Summary Statistics")
    print("=" * 60)
    
    print("\nDrought Events Summary:")
    print(f"  Total events: {len(result_compound['drought_events'])}")
    print(f"  Mean duration: {result_compound['drought_events']['duration'].mean():.1f} days")
    print(f"  Mean severity: {result_compound['drought_events']['severity'].mean():.2f}")
    print(f"  Max severity: {result_compound['drought_events']['severity'].max():.2f}")
    
    print("\nHeatwave Events Summary:")
    print(f"  Total events: {len(result_compound['heatwave_events'])}")
    print(f"  Mean duration: {result_compound['heatwave_events']['duration'].mean():.1f} days")
    print(f"  Mean intensity: {result_compound['heatwave_events']['intensity'].mean():.2f}")
    print(f"  Max intensity: {result_compound['heatwave_events']['intensity'].max():.2f}")
    
    print("\nCompound Events Summary:")
    print(f"  Total compound days: {result_compound['compound_daily']['is_compound'].sum()}")
    n_compound_events = result_compound['compound_flags']['compound_order'].max()
    print(f"  Number of compound events: {int(n_compound_events)}")
    
    print("\n" + "=" * 60)
    print("Demo completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
