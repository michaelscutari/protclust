# Contributing to protclust

Thank you for considering contributing to protclust! This document provides guidelines and instructions for contributing.

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/username/protclust.git
   cd protclust
   ```

2. Create a virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -e .
   pip install pytest pytest-cov pre-commit ruff
   ```

3. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

## Running Tests

Run tests using pytest:

```bash
pytest
```

For test coverage information:

```bash
pytest --cov=protclust tests/
```

## Code Style

This project uses `ruff` for code formatting and linting. Pre-commit hooks will automatically check your code before commits.

To manually run ruff:

```bash
ruff check .
ruff format .
```

## Pull Request Process

1. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes and commit them with descriptive commit messages.

3. Ensure all tests pass and code style checks pass.

4. Push your branch to GitHub:
   ```bash
   git push origin feature/your-feature-name
   ```

5. Open a pull request against the `main` branch and provide a clear description of the changes.

## Reporting Issues

If you find a bug or have a feature request, please create an issue on the GitHub repository. Include as much detail as possible:

- For bugs: steps to reproduce, expected behavior, and actual behavior
- For features: clear description of the proposed functionality and use cases

## Adding New Features

When adding new features, please:

1. Include appropriate tests
2. Update documentation
3. Consider backwards compatibility
4. Update type hints where applicable

Thank you for your contributions!
