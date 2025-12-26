"""Tests for leap year normalization (to_365)."""

import numpy as np
import pandas as pd
import pytest

from Shan_daily import to_365


def test_to_365_basic():
    """Test basic leap year normalization."""
    # Create data with a leap year (2000)
    dates = pd.date_range("2000-01-01", "2000-12-31", freq="D")
    data = pd.DataFrame({
        "date": dates,
        "value": np.arange(len(dates))
    })
    
    result = to_365(data, how="mean")
    
    # Should have 365 days after normalization
    assert len(result) == 365
    assert "value" in result.columns


def test_to_365_sum_method():
    """Test to_365 with sum aggregation."""
    dates = pd.date_range("2000-02-28", "2000-03-02", freq="D")
    data = pd.DataFrame({
        "date": dates,
        "pre": [10, 20, 30, 40]  # Feb 28, 29, Mar 1, 2
    })
    
    result = to_365(data, how="sum")
    
    # Feb 29 and Mar 1 should be combined
    assert len(result) == 3


def test_to_365_mean_method():
    """Test to_365 with mean aggregation."""
    dates = pd.date_range("2000-02-28", "2000-03-02", freq="D")
    data = pd.DataFrame({
        "date": dates,
        "tem": [10, 20, 30, 40]
    })
    
    result = to_365(data, how="mean")
    
    # Feb 29 and Mar 1 should be averaged
    assert len(result) == 3


def test_to_365_no_leap_year():
    """Test that non-leap years pass through unchanged."""
    dates = pd.date_range("2001-01-01", "2001-12-31", freq="D")
    data = pd.DataFrame({
        "date": dates,
        "value": np.arange(365)
    })
    
    result = to_365(data, how="mean")
    
    assert len(result) == 365
    np.testing.assert_array_almost_equal(result["value"], data["value"])

