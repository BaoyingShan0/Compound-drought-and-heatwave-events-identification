"""Pytest configuration and shared fixtures."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_daily_data():
    """Generate sample daily climate data for testing."""
    dates = pd.date_range("2000-01-01", periods=365*5, freq="D")
    
    # Generate synthetic precipitation data (0-50 mm)
    np.random.seed(42)
    pre = np.random.gamma(2, 3, size=len(dates))
    pre[pre < 0] = 0
    
    pre_df = pd.DataFrame({
        "date": dates,
        "pre": pre
    })
    
    # Generate synthetic temperature data (10-30°C with seasonal cycle)
    day_of_year = dates.dayofyear
    seasonal = 10 * np.sin(2 * np.pi * (day_of_year - 80) / 365)
    noise = np.random.normal(0, 2, size=len(dates))
    tem_mean = 20 + seasonal + noise
    
    tem_df = pd.DataFrame({
        "date": dates,
        "tem_mean": tem_mean
    })
    
    return pre_df, tem_df


@pytest.fixture
def sample_standardized_index():
    """Generate sample standardized index values."""
    np.random.seed(42)
    dates = pd.date_range("2000-01-01", periods=365*5, freq="D")
    si_values = np.random.normal(0, 1, size=len(dates))
    
    return pd.Series(si_values, index=dates)

