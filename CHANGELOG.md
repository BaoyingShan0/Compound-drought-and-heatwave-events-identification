# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Project structure improvements
- Standard documentation files

## [0.1.0] - 2024-12-26

### Added
- Initial Python implementation of compound drought-heatwave event identification
- Main entry functions: `identify_compound_events` and `identify_extremes`
- Standardized indices computation (SI_nonparametric and SI_best_distribution)
- PRM (Pre-identification, Removal, Merging) algorithm for event identification
- Threshold optimization via grid search with statistical validation
- Support for four compound event types (AND, OR, conditional combinations)
- Leap year normalization functionality
- Visualization tools for extreme event identification
- Comprehensive example scripts and demonstration data
- Support for pandas and xarray data structures
- Mann-Kendall trend test for non-stationarity detection

### Documentation
- Comprehensive README with usage examples
- API documentation in code
- Example demonstrations in `examples/run_demo.py`
- MATLAB version preserved in `Matlab_version/` for reference

### Dependencies
- Core: numpy, pandas, scipy, matplotlib
- Optional: numba, xarray, pymannkendall, bottleneck

[Unreleased]: https://github.com/BaoyingShan0/Compound-drought-and-heatwave-events-identification/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/BaoyingShan0/Compound-drought-and-heatwave-events-identification/releases/tag/v0.1.0
