<p align="left">
  <img src="assets/images/logo.png" alt="protclust logo" width="100"/>
</p>

# protclust

[![PyPI version](https://img.shields.io/pypi/v/protclust.svg)](https://pypi.org/project/protclust/)
[![Tests](https://github.com/michaelscutari/protclust/workflows/Tests/badge.svg)](https://github.com/michaelscutari/protclust/actions)
[![Coverage](https://img.shields.io/badge/Coverage-85%25-green)](https://github.com/YOUR-USERNAME/protclust/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/pypi/pyversions/protclust.svg)](https://pypi.org/project/protclust/)

A Python library for working with protein sequence data, providing:
- Clustering capabilities via MMseqs2
- Machine learning dataset creation with cluster-aware splits

---

## Requirements

This library requires [MMseqs2](https://github.com/soedinglab/MMseqs2), which must be installed and accessible via the command line. MMseqs2 can be installed using one of the following methods:

### Installation Options for MMseqs2

- **Homebrew**:
    ```bash
    brew install mmseqs2
    ```

- **Conda**:
    ```bash
    conda install -c conda-forge -c bioconda mmseqs2
    ```

- **Docker**:
    ```bash
    docker pull ghcr.io/soedinglab/mmseqs2
    ```

- **Static Build (AVX2, SSE4.1, or SSE2)**:
    ```bash
    wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
    tar xvfz mmseqs-linux-avx2.tar.gz
    export PATH=$(pwd)/mmseqs/bin/:$PATH
    ```

MMseqs2 must be accessible via the `mmseqs` command in your system's PATH. If the library cannot detect MMseqs2, it will raise an error.

## Installation

### Installation

You can install protclust using pip:

```bash
pip install protclust
```

Or if installing from source, clone the repository and run:

```bash
pip install -e .
```

For development purposes, also install the testing dependencies:

```bash
pip install pytest pytest-cov pre-commit ruff
```

## Features

### Sequence Clustering and Dataset Creation

```python
import pandas as pd
from protclust import clean, cluster, split, set_verbosity

# Enable detailed logging (optional)
set_verbosity(verbose=True)

# Example data
df = pd.DataFrame({
    "id": ["seq1", "seq2", "seq3", "seq4"],
    "sequence": ["ACDEFGHIKL", "ACDEFGHIKL", "MNPQRSTVWY", "MNPQRSTVWY"]
})

# Clean data
clean_df = clean(df, sequence_col="sequence")

# Cluster sequences
clustered_df = cluster(clean_df, sequence_col="sequence", id_col="id")

# Split data into train and test sets
train_df, test_df = split(clustered_df, group_col="cluster_representative", test_size=0.3)

print("Train set:\n", train_df)
print("Test set:\n", test_df)

# MILP-based splitting with property balancing
from protclust import milp_split
train_df, test_df = milp_split(
    clustered_df,
    group_col="cluster_representative",
    test_size=0.3,
    balance_cols=["molecular_weight", "hydrophobicity"]
)
```

## Parameters

Common parameters for clustering functions:

- `df`: Pandas DataFrame containing sequence data
- `sequence_col`: Column name containing sequences
- `id_col`: Column name containing unique identifiers
- `min_seq_id`: Minimum sequence identity threshold (0.0-1.0, default 0.3)
- `coverage`: Minimum alignment coverage (0.0-1.0, default 0.5)
- `cov_mode`: Coverage mode (0-3, default 0)
- `cluster_mode`: Clustering algorithm (0: Set-Cover, 1: Connected component, 2: Greedy by length, default 0)
- `cluster_steps`: Number of cascaded clustering steps for large datasets (default 1)
- `test_size`: Desired fraction of data in test set (default 0.2)
- `random_state`: Random seed for reproducibility
- `tolerance`: Acceptable deviation from desired split sizes (default 0.05)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Run tests (`pytest tests/`)
4. Commit your changes (`git commit -m 'Add some amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use protclust in your research, please cite:

```bibtex
@software{protclust,
  author = {Michael Scutari},
  title = {protclust: Protein Sequence Clustering and ML Dataset Creation},
  url = {https://github.com/michaelscutari/protclust},
  version = {0.2.0},
  year = {2025},
}
```
