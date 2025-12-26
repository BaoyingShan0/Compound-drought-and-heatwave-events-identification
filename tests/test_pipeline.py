"""Tests for main pipeline functions."""

import pandas as pd
import pytest

from Shan_daily import identify_compound_events, identify_extremes


def test_identify_extremes(sample_daily_data):
    """Test single-variable extreme identification."""
    pre_df, _ = sample_daily_data
    
    result = identify_extremes(
        climate_daily=pre_df,
        extreme_type="d",
        scale=30,
        group_by="sum",
        si_method="nonparametric",
        nsp=None
    )
    
    # Check output structure
    assert isinstance(result, dict)
    assert "extreme_events" in result
    assert "extreme_daily" in result
    assert "standardized_index" in result
    
    # Check events dataframe
    events = result["extreme_events"]
    assert isinstance(events, pd.DataFrame)
    assert "duration" in events.columns
    assert "severity" in events.columns
    assert "intensity" in events.columns


def test_identify_compound_events(sample_daily_data):
    """Test compound event identification."""
    pre_df, tem_df = sample_daily_data
    
    result = identify_compound_events(
        pre_daily=pre_df,
        tem_daily=tem_df,
        scale_pre=30,
        scale_tem=3,
        nsp=None,
        si_method="nonparametric",
        compound_type=1
    )
    
    # Check output structure
    assert isinstance(result, dict)
    assert "drought_events" in result
    assert "heatwave_events" in result
    assert "compound_daily" in result
    assert "spi" in result
    assert "shi" in result
    
    # Check compound daily dataframe
    compound_daily = result["compound_daily"]
    assert "is_compound" in compound_daily.columns
    assert compound_daily["is_compound"].dtype in [int, bool]


def test_identify_extremes_heatwave(sample_daily_data):
    """Test heatwave identification."""
    _, tem_df = sample_daily_data
    
    result = identify_extremes(
        climate_daily=tem_df,
        extreme_type="h",
        scale=3,
        group_by="mean",
        si_method="nonparametric"
    )
    
    assert isinstance(result, dict)
    assert "extreme_events" in result
    events = result["extreme_events"]
    assert isinstance(events, pd.DataFrame)


def test_compound_types(sample_daily_data):
    """Test different compound event types."""
    pre_df, tem_df = sample_daily_data
    
    for compound_type in [1, 2, 3, 4]:
        result = identify_compound_events(
            pre_daily=pre_df,
            tem_daily=tem_df,
            scale_pre=30,
            scale_tem=3,
            compound_type=compound_type,
            si_method="nonparametric"
        )
        
        assert isinstance(result, dict)
        assert "compound_daily" in result

