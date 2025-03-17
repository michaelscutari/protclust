# Changelog

All notable changes to protclust will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-03-17

### Added
- Initial release of protclust with the following features:
  - Sequence preprocessing and cleaning
  - MMseqs2 integration for protein sequence clustering
  - Cluster-aware train/test splitting
  - Multiway (train/val/test) splitting
  - Constrained splitting with custom ID assignments
  - K-fold cross-validation with cluster awareness
  - MILP-based splitting with property balancing
  - Embedding support for proteins:
    - Baseline embedders (BLOSUM62, AAC, Property, One-hot)
    - Integration with ESM models
    - ProtTrans embeddings
    - Remote API embeddings
  - Embedding utility functions (dimension reduction, storage)
