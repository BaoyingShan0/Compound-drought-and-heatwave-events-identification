# Examples

This directory contains example scripts and data for demonstrating the Shan_daily package functionality.

## 📁 Contents

- **`run_demo.py`**: Comprehensive demonstration script covering all major features
- **`data/`**: Example climate data files
  - `daily_pre.csv`: Daily precipitation data (mm)
  - `tem_daily.csv`: Daily mean temperature data (°C)

## 🚀 Quick Start

### Run the Complete Demo

```bash
python run_demo.py
```

This script demonstrates:

1. **Complete Workflow (Main Entry Functions)**
   - Compound drought-heatwave event identification
   - Single-variable extreme event identification

2. **Individual Component Usage**
   - Data preprocessing (leap year normalization)
   - Standardized indices computation (SPI/SHI)
   - Threshold optimization
   - Event identification with PRM algorithm
   - Compound event detection

3. **Visualization**
   - Plotting extreme event identification results

4. **Different Compound Types**
   - Type 1: AND (intersection)
   - Type 2: OR (union)
   - Type 3: Drought-conditional-heatwave
   - Type 4: Heatwave-conditional-drought

## 📊 Example Data

### Precipitation Data Format (`daily_pre.csv`)

| date       | pre  |
|------------|------|
| 1950-01-01 | 2.5  |
| 1950-01-02 | 0.0  |
| 1950-01-03 | 5.3  |
| ...        | ...  |

### Temperature Data Format (`tem_daily.csv`)

| date       | tem_mean |
|------------|----------|
| 1950-01-01 | 15.2     |
| 1950-01-02 | 16.8     |
| 1950-01-03 | 14.5     |
| ...        | ...      |

## 💡 Usage Tips

### Customize Parameters

Edit `run_demo.py` to experiment with different parameters:

```python
# Change time scales
result = identify_compound_events(
    pre_daily=pre_daily,
    tem_daily=tem_daily,
    scale_pre=60,    # Try 60-day SPI instead of 30
    scale_tem=7,     # Try 7-day SHI instead of 3
    ...
)
```

### Use Your Own Data

Replace the data files with your own:

1. Prepare CSV files with the same format
2. Update file paths in the script:

```python
pre_raw = pd.read_csv("path/to/your/daily_pre.csv")
tem_raw = pd.read_csv("path/to/your/tem_daily.csv")
```

### Save Results

Uncomment the save lines in `run_demo.py` to export results:

```python
# Save results
drought_events.to_csv("outputs/drought_events.csv", index=False)
heatwave_events.to_csv("outputs/heatwave_events.csv", index=False)
compound_daily.to_csv("outputs/compound_daily.csv", index=False)
```

## 📈 Expected Output

The demo will print:

- Data dimensions and statistics
- Number of identified events
- Event characteristics (duration, severity, intensity)
- Statistical validation results
- Optimal thresholds
- Compound event counts

Visualization files (if saved):
- `drought_identification.png`
- `heatwave_identification.png`

## 🔧 Troubleshooting

### ImportError

If you get import errors, make sure the package is installed:

```bash
cd ..
pip install -e .
```

### Memory Issues

For large datasets, consider:
- Using the non-parametric method (`si_method="nonparametric"`)
- Processing data in chunks
- Reducing the window size (`nsp` parameter)

### Performance Optimization

To speed up computation:
- Install optional dependencies: `pip install numba pymannkendall xarray`
- Use `si_method="nonparametric"` (faster than `"best_distribution"`)
- Reduce grid search resolution in threshold optimization

## 📖 More Information

For detailed API documentation and methodology, see the main [README.md](../README.md) in the project root.

## 🤝 Contributing

Found an issue or have an example to share? See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

