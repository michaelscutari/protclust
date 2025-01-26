# My Library

A Python library for clustering biological sequences and splitting data into train/test sets. This library uses MMseqs2 for sequence clustering.

---

## Requirements

This library requires [MMseqs2](https://github.com/soedinglab/MMseqs2), which must be installed and accessible via the command line. MMseqs2 can be installed using one of the following methods:

### Installation Options

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

MMseqs2 must be accessible via the `mmseqs` command in your system’s PATH. If the library cannot detect MMseqs2, it will raise an error.

## Installation

To install this library, clone the repository and run the following command in the project root:

```bash
pip install .
```

## Features

This library provides two main functions:

1. **cluster_sequences**: Clusters biological sequences using MMseqs2 and adds a `representative_sequence` column to a pandas DataFrame.
2. **cluster_split_data**: Clusters biological sequences and splits the data into train/test sets. The split ensures no overlap between clusters in the training and testing sets, maintaining cluster integrity.

## Quick Start

Here’s an example of how to use this library:

```python
import pandas as pd
from my_library import cluster_sequences, cluster_split_data

# Example data
df = pd.DataFrame({
        "id": ["seq1", "seq2", "seq3"],
        "sequence": ["ATCG", "GCTA", "TAGC"]
})

# Clean data
clean_df = clean(df)

# Cluster sequences
clustered_df = cluster(clean_df, sequence_col="sequence")

# Split data into train and test sets
train_df, test_df = cluster_split(clustered_df, sequence_col="sequence", test_size=0.3)

print("Train set:\n", train_df)
print("Test set:\n", test_df)
```

## Notes

- Ensure MMseqs2 is installed and accessible via the command line before using the library.
- You can customize MMseqs2 parameters such as `min_seq_id`, `coverage`, and `cov_mode` when calling the functions.