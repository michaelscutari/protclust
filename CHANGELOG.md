# Changelog

All notable changes to protclust will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.5] - 2025-03-21

### Enhancements
- Added Raygun embeddings.

### Fixed
- Some tests still called `cluster` using the removed `threads` parameter.

## [0.1.4] - 2025-03-20
### Enhancements

- **Reproducible Clustering**: Added deterministic behavior to MMseqs2 clustering

New random_state parameter for consistent results across runs
Added threads parameter (default: 1) for deterministic execution
Disabled sequence shuffling in MMseqs2 (--shuffle 0)
Added input sequence sorting for consistent processing order

### Fixed

- Eliminated non-deterministic behavior in sequence clustering

## [0.1.2] - 2025-03-19

### Enhancements

- **Improved Logging System**:
 - Added a more flexible `set_verbosity()` function that accepts numeric levels (0=WARNING, 1=INFO, 2=DEBUG) for finer control over log output
 - Added `get_verbosity()` function to retrieve the current logging level
 - Changed default logging level to WARNING to reduce unnecessary console output in production environments

- **MILP Solver Silence Control**: The MILP solver now respects logging levels, only showing detailed output at INFO level or above

### Fixed
- Fixed logging configuration to properly initialize at startup
- Improved console output clarity by controlling verbosity of external libraries

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
