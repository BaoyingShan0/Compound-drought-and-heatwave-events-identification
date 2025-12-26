"""Tests for non-parametric standardized index computation."""

import numpy as np
import pandas as pd
import pytest

from Shan_daily import SI_nonparametric


def test_SI_nonparametric_basic(sample_daily_data):
    """Test basic SI_nonparametric computation."""
    pre_df, _ = sample_daily_data
    
    result = SI_nonparametric(
        pre_df["pre"],
        scale=30,
        group_by="sum",
        NSP=None,
        pp_method="gringorten"
    )
    
    # Result should be a Series
    assert isinstance(result, pd.Series)
    
    # Should have values in reasonable range for standardized index
    assert result.min() > -5
    assert result.max() < 5


def test_SI_nonparametric_stationary():
    """Test SI computation with stationary assumption."""
    np.random.seed(42)
    dates = pd.date_range("2000-01-01", periods=365*10, freq="D")
    values = np.random.gamma(2, 3, size=len(dates))
    
    result = SI_nonparametric(
        pd.Series(values, index=dates),
        scale=1,
        group_by="sum",
        NSP=None
    )
    
    # Mean should be close to 0
    assert abs(result.mean()) < 0.2
    
    # Std should be close to 1
    assert abs(result.std() - 1.0) < 0.3


def test_SI_nonparametric_nonstationary():
    """Test SI computation with non-stationary window."""
    np.random.seed(42)
    dates = pd.date_range("2000-01-01", periods=365*30, freq="D")
    values = np.random.gamma(2, 3, size=len(dates))
    
    result = SI_nonparametric(
        pd.Series(values, index=dates),
        scale=30,
        group_by="sum",
        NSP=30  # 30-year window
    )
    
    assert isinstance(result, pd.Series)
    assert len(result) == len(values)

