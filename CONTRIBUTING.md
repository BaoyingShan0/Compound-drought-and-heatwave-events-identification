# Contributing to Compound Drought and Heatwave Events Identification

Thank you for your interest in contributing to this project! 🎉

## 📋 Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

This project adheres to a code of conduct that all contributors are expected to follow. Please be respectful and constructive in all interactions.

## How Can I Contribute?

### Reporting Bugs

If you find a bug, please open an issue with:
- A clear, descriptive title
- Steps to reproduce the problem
- Expected vs. actual behavior
- Your environment (Python version, OS, package versions)
- Minimal code example demonstrating the issue

### Suggesting Enhancements

Enhancement suggestions are welcome! Please open an issue describing:
- The motivation for the enhancement
- How it would be used
- Potential implementation approach (if you have ideas)

### Contributing Code

1. **Fork the repository** and create a new branch from `main`
2. **Make your changes** following the coding standards below
3. **Add tests** if applicable
4. **Update documentation** if you're changing functionality
5. **Submit a pull request** with a clear description of your changes

## Development Setup

### 1. Clone the Repository

```bash
git clone https://github.com/BaoyingShan0/Compound-drought-and-heatwave-events-identification.git
cd Compound-drought-and-heatwave-events-identification
```

### 2. Create a Development Environment

Using conda (recommended):

```bash
conda env create -f environment.yml
conda activate shan_daily
```

Or using pip:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e ".[full]"
```

### 3. Install Development Dependencies

```bash
pip install pytest black flake8 mypy
```

## Coding Standards

### Python Style

- Follow [PEP 8](https://pep8.org/) style guide
- Use **4 spaces** for indentation (no tabs)
- Maximum line length: **88 characters** (Black default)
- Use descriptive variable names

### Code Formatting

Format your code with Black before committing:

```bash
black src/Shan_daily/
```

### Type Hints

Add type hints to function signatures when possible:

```python
def process_data(values: np.ndarray, scale: int = 30) -> pd.DataFrame:
    ...
```

### Documentation

- Add docstrings to all public functions and classes
- Use NumPy-style docstrings:

```python
def example_function(param1: int, param2: str) -> bool:
    """
    Brief description of the function.
    
    Parameters
    ----------
    param1 : int
        Description of param1
    param2 : str
        Description of param2
        
    Returns
    -------
    bool
        Description of return value
        
    Examples
    --------
    >>> example_function(5, "test")
    True
    """
    ...
```

### Testing

- Write tests for new functionality
- Place tests in the `tests/` directory
- Use descriptive test names: `test_<function>_<scenario>`
- Run tests before submitting:

```bash
pytest tests/
```

## Submitting Changes

### Commit Messages

Write clear, concise commit messages:

```
Add feature for multi-scale SI computation

- Implement rolling window aggregation
- Add support for custom time scales
- Update documentation with examples
```

Format:
- First line: brief summary (50 chars or less)
- Blank line
- Detailed description with bullet points if needed

### Pull Request Process

1. **Update documentation** if you're changing functionality
2. **Update CHANGELOG.md** under the "Unreleased" section
3. **Ensure tests pass** and code is formatted
4. **Create a pull request** with:
   - Clear title describing the change
   - Description of what and why
   - Link to related issues (if any)
5. **Respond to feedback** from reviewers

### PR Review Criteria

Your PR will be reviewed for:
- ✅ Code quality and style
- ✅ Test coverage
- ✅ Documentation completeness
- ✅ Backward compatibility
- ✅ Performance considerations

## Development Tips

### Running Examples

Test your changes with the demo:

```bash
cd examples
python run_demo.py
```

### Checking Code Quality

```bash
# Format code
black src/Shan_daily/

# Check style
flake8 src/Shan_daily/

# Type checking
mypy src/Shan_daily/
```

### Building Documentation

If you're updating documentation:

```bash
cd docs
make html
```

## Questions?

If you have questions about contributing, feel free to:
- Open an issue with your question
- Contact the maintainer: baoying.shan@polimi.it

## Recognition

All contributors will be acknowledged in the project documentation. Thank you for helping improve this project! 🙌

