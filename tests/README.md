# Tests

This directory contains the test suite for the Shan_daily package.

## 🧪 Running Tests

### Run all tests

```bash
pytest
```

### Run with coverage

```bash
pytest --cov=Shan_daily --cov-report=html
```

### Run specific test file

```bash
pytest tests/test_pipeline.py
```

### Run specific test function

```bash
pytest tests/test_pipeline.py::test_identify_extremes
```

### Run with verbose output

```bash
pytest -v
```

## 📁 Test Structure

- `conftest.py`: Shared fixtures and pytest configuration
- `test_to_365.py`: Tests for leap year normalization
- `test_SI_nonparametric.py`: Tests for standardized index computation
- `test_pipeline.py`: Tests for main entry functions
- `README.md`: This file

## 🔧 Test Fixtures

Common fixtures are defined in `conftest.py`:

- `sample_daily_data`: Synthetic daily precipitation and temperature data
- `sample_standardized_index`: Synthetic standardized index values

## ✅ Test Coverage Goals

Target test coverage areas:

- ✅ Data preprocessing (leap year normalization)
- ✅ Standardized index computation (parametric and non-parametric)
- ✅ Event identification (PRM algorithm)
- ✅ Threshold optimization
- ✅ Compound event detection
- ✅ Main pipeline functions

## 📝 Writing New Tests

When adding new tests:

1. Follow the naming convention: `test_<module>_<functionality>`
2. Use descriptive test function names: `test_function_with_specific_scenario`
3. Add docstrings to explain what is being tested
4. Use fixtures from `conftest.py` when possible
5. Test both happy paths and edge cases

Example:

```python
def test_new_functionality():
    """Test that new_functionality works correctly."""
    # Arrange
    input_data = ...
    
    # Act
    result = new_functionality(input_data)
    
    # Assert
    assert result.shape == expected_shape
    assert result.columns == expected_columns
```

## 🐛 Testing Tips

- Use `pytest.mark.parametrize` for testing multiple inputs
- Use `pytest.raises` for testing expected exceptions
- Keep tests isolated and independent
- Mock external dependencies when needed
- Test edge cases (empty data, NaN values, etc.)

## 📊 Continuous Integration

Tests should pass before merging pull requests. Configure CI to run tests automatically on push.

