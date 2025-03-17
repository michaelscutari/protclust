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
- Protein sequence embeddings and feature extraction

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
train_df, test_df = split(clustered_df, group_col="representative_sequence", test_size=0.3)

print("Train set:\n", train_df)
print("Test set:\n", test_df)

# Or use the combined function
from protclust import train_test_cluster_split
train_df, test_df = train_test_cluster_split(df, sequence_col="sequence", id_col="id", test_size=0.3)
```

### Advanced Splitting Options

```python
# Three-way split
from protclust import train_test_val_cluster_split
train_df, val_df, test_df = train_test_val_cluster_split(
    df, sequence_col="sequence", test_size=0.2, val_size=0.1
)

# K-fold cross-validation with cluster awareness
from protclust import cluster_kfold
folds = cluster_kfold(df, sequence_col="sequence", n_splits=5)

# MILP-based splitting with property balancing
from protclust import milp_split
train_df, test_df = milp_split(
    clustered_df,
    group_col="representative_sequence",
    test_size=0.3,
    balance_cols=["molecular_weight", "hydrophobicity"]
)
```

### Protein Embeddings

```python
# Basic embeddings
from protclust.embeddings import blosum62, aac, property_embedding, onehot

# Add BLOSUM62 embeddings
df_with_blosum = blosum62(df, sequence_col="sequence")

# Generate embeddings with ESM models (requires extra dependencies)
from protclust.embeddings import embed_sequences

# ESM embedding
df_with_esm = embed_sequences(df, "esm", sequence_col="sequence")

# Saving embeddings to HDF5
df_with_refs = embed_sequences(
    df,
    "esm",
    sequence_col="sequence",
    use_hdf=True,
    hdf_path="embeddings.h5"
)

# Retrieve embeddings
from protclust.embeddings import get_embeddings
embeddings = get_embeddings(df_with_refs, "esm", hdf_path="embeddings.h5")
```

## Parameters

Common parameters for clustering functions:

- `df`: Pandas DataFrame containing sequence data
- `sequence_col`: Column name containing sequences
- `id_col`: Column name containing unique identifiers
- `min_seq_id`: Minimum sequence identity threshold (0.0-1.0, default 0.3)
- `coverage`: Minimum alignment coverage (0.0-1.0, default 0.5)
- `cov_mode`: Coverage mode (0-3, default 0)
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
  version = {0.1.0},
  year = {2025},
}
```
