import logging
import os
import random
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

# Set up logging
logger = logging.getLogger("protein_clustering")
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)  # Default level

# default seed
_DEFAULT_SEED = 42
_current_seed = _DEFAULT_SEED


def set_seed(seed=_DEFAULT_SEED):
    """
    Set random seeds for reproducibility.

    Parameters:
        seed (int): Random seed value (default 42)
    """
    global _current_seed
    if seed != _current_seed:
        logger.info(f"Setting random seed to {seed}")
        random.seed(seed)
        np.random.seed(seed)
        _current_seed = seed


def set_verbosity(verbose=False):
    """
    Set the verbosity level for the package.

    Parameters:
        verbose (bool or int): If True or 1, sets to INFO level.
                               If 2, sets to DEBUG level.
                               If False or 0, sets to WARNING level.
    """
    if verbose == 2:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbosity set to DEBUG level")
    elif verbose:
        logger.setLevel(logging.INFO)
        logger.info("Verbosity set to INFO level")
    else:
        logger.setLevel(logging.WARNING)


def _check_mmseqs():
    """
    Ensures 'mmseqs' command is in PATH.
    """
    logger.debug("Checking if MMseqs2 is installed")
    if shutil.which("mmseqs") is None:
        logger.error("MMseqs2 not found in PATH")
        raise EnvironmentError(
            "MMseqs2 is not installed or not found in PATH. "
            "See the README for installation instructions."
        )
    logger.debug("MMseqs2 found in PATH")


def clean(df, sequence_col="sequence", valid_amino_acids="ACDEFGHIKLMNPQRSTVWY"):
    """
    Removes sequences with invalid protein characters.

    Parameters:
        df (pd.DataFrame): Input DataFrame with protein sequences.
        sequence_col (str): Name of the column containing sequences.
        valid_amino_acids (str): String of valid amino acid characters.

    Returns:
        pd.DataFrame: Cleaned DataFrame with only valid sequences.
    """
    logger.info(
        f"Cleaning sequences in column '{sequence_col}' with valid amino acids: {valid_amino_acids}"
    )
    logger.info(f"Input dataframe has {len(df)} sequences")

    df[sequence_col] = df[sequence_col].str.upper()
    df = df.dropna(subset=[sequence_col])

    logger.debug(f"After removing NaN values: {len(df)} sequences")

    valid_sequence_mask = df[sequence_col].apply(
        lambda seq: all(aa in valid_amino_acids for aa in seq)
    )

    result_df = df[valid_sequence_mask].reset_index(drop=True)

    invalid_count = len(df) - len(result_df)
    logger.info(f"Removed {invalid_count} sequences with invalid amino acids")
    logger.info(f"Final dataframe has {len(result_df)} valid sequences")

    return result_df


def cluster(
    df,
    sequence_col,
    id_col=None,
    min_seq_id=0.3,
    coverage=0.5,
    cov_mode=0,
    alignment_mode=0,
):
    """
    Clusters sequences with MMseqs2 and adds a 'representative_sequence' column.

    Parameters:
        df (pd.DataFrame): Input DataFrame with columns for IDs and sequences.
        sequence_col (str): Name of the column containing sequences.
        id_col (str): Unique ID column (default "id").
        min_seq_id (float): Minimum sequence identity for clustering (default 0.3).
        coverage (float): Minimum alignment coverage (default 0.5).
        cov_mode (int): Coverage mode for MMseqs2 (default 0).
        alignment_mode (int): Alignment mode for MMseqs2 (default 0).

    Returns:
        pd.DataFrame: Original DataFrame with a new 'representative_sequence' column.
    """
    logger.info("Starting sequence clustering with MMseqs2")
    logger.info(
        f"Parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}, alignment_mode={alignment_mode}"
    )

    _check_mmseqs()

    if id_col is None:
        df = df.reset_index()
        id_col = "index"
        logger.debug(f"No id_col provided, using '{id_col}' as identifier")

    if sequence_col not in df or id_col not in df:
        logger.error(f"Required columns missing: {sequence_col} or {id_col}")
        raise ValueError(f"The DataFrame must have '{id_col}' and '{sequence_col}'.")

    logger.info(f"Clustering {len(df)} sequences")

    df["sanitized_id"] = df[id_col].str.replace(" ", "_")
    tmp_dir = tempfile.mkdtemp()
    logger.debug(f"Created temporary directory: {tmp_dir}")

    try:
        input_fasta = os.path.join(tmp_dir, "input.fasta")
        with open(input_fasta, "w") as fasta_file:
            for _, row in df.iterrows():
                fasta_file.write(f">{row['sanitized_id']}\n{row[sequence_col]}\n")

        logger.debug(f"Wrote {len(df)} sequences to FASTA file")

        output_dir = os.path.join(tmp_dir, "output")
        tmp_mmseqs = os.path.join(tmp_dir, "tmp_mmseqs")

        mmseqs_cmd = [
            "mmseqs",
            "easy-cluster",
            input_fasta,
            output_dir,
            tmp_mmseqs,
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
            "--alignment-mode",
            str(alignment_mode),
        ]

        logger.debug(f"Running MMseqs2 command: {' '.join(mmseqs_cmd)}")

        if logger.level <= logging.DEBUG:
            subprocess.run(mmseqs_cmd, check=True)
        else:
            subprocess.run(
                mmseqs_cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        clusters_file = os.path.join(output_dir + "_cluster.tsv")
        if not os.path.exists(clusters_file):
            logger.error("MMseqs2 clustering results file not found")
            raise FileNotFoundError("MMseqs2 clustering results not found.")

        logger.debug(f"Reading clustering results from {clusters_file}")

        cluster_map = {}
        cluster_sizes = {}
        with open(clusters_file, "r") as f:
            for line in f:
                rep, seq = line.strip().split("\t")
                cluster_map[seq] = rep
                cluster_sizes[rep] = cluster_sizes.get(rep, 0) + 1

        logger.info(f"Found {len(cluster_sizes)} clusters")

        if logger.level <= logging.DEBUG:
            # Report cluster distribution statistics
            cluster_size_counts = {}
            for size in cluster_sizes.values():
                cluster_size_counts[size] = cluster_size_counts.get(size, 0) + 1

            logger.debug("Cluster size distribution:")
            for size in sorted(cluster_size_counts.keys()):
                logger.debug(f"  Size {size}: {cluster_size_counts[size]} clusters")

        reverse_map = dict(zip(df["sanitized_id"], df[id_col]))
        df["representative_sequence"] = df["sanitized_id"].apply(
            lambda x: reverse_map.get(cluster_map.get(x, x), x)
        )

        logger.info(
            "Clustering complete, added 'representative_sequence' column to DataFrame"
        )

    finally:
        logger.debug(f"Cleaning up temporary directory: {tmp_dir}")
        shutil.rmtree(tmp_dir)

    df.drop(columns=["sanitized_id"], inplace=True)
    return df


def split(
    df,
    group_col="representative_sequence",
    test_size=0.2,
    random_state=None,
    tolerance=0.05,
):
    """
    Splits DataFrame into train/test sets based on grouping in a specified column.

    Parameters:
        df (pd.DataFrame): DataFrame to split.
        group_col (str): Column by which to group before splitting.
        test_size (float): Desired fraction of data in test set (default 0.2).
        random_state (int): Optional random state for reproducibility (unused in subset-sum).
        tolerance (float): Acceptable deviation from test_size (default 0.05).

    Returns:
        (pd.DataFrame, pd.DataFrame): (train_df, test_df)
    """
    logger.info(f"Splitting data by '{group_col}' with target test size {test_size}")

    total_sequences = len(df)
    target_test_count = int(round(test_size * total_sequences))

    logger.info(f"Total sequence count: {total_sequences}")
    logger.info(f"Target test count: {target_test_count}")

    group_sizes_df = df.groupby(group_col).size().reset_index(name="group_size")

    logger.debug(f"Found {len(group_sizes_df)} unique groups in '{group_col}'")

    groups = group_sizes_df[group_col].tolist()
    sizes = group_sizes_df["group_size"].tolist()

    logger.debug("Finding optimal subset-sum solution for test set")

    dp = {0: []}
    for idx, group_size in enumerate(sizes):
        current_dp = dict(dp)
        for current_sum, idx_list in dp.items():
            new_sum = current_sum + group_size
            if new_sum not in current_dp:
                current_dp[new_sum] = idx_list + [idx]
        dp = current_dp

    best_sum = min(dp.keys(), key=lambda s: abs(s - target_test_count))
    best_group_indices = dp[best_sum]
    chosen_groups = [groups[i] for i in best_group_indices]

    logger.debug(f"Best achievable test set size: {best_sum} sequences")
    logger.debug(f"Selected {len(chosen_groups)} groups for test set")

    test_df = df[df[group_col].isin(chosen_groups)]
    train_df = df[~df[group_col].isin(chosen_groups)]

    achieved_test_fraction = len(test_df) / total_sequences

    logger.info(
        f"Train set: {len(train_df)} sequences ({len(train_df) / total_sequences:.2%})"
    )
    logger.info(f"Test set: {len(test_df)} sequences ({achieved_test_fraction:.2%})")

    if abs(achieved_test_fraction - test_size) > tolerance:
        logger.warning(
            f"Desired test fraction = {test_size:.2f}, "
            f"achieved = {achieved_test_fraction:.2f}. "
            "This is the closest possible split given the constraint to keep groups together."
        )

    return train_df, test_df


def train_test_cluster_split(
    df,
    sequence_col,
    id_col=None,
    test_size=0.2,
    min_seq_id=0.3,
    coverage=0.5,
    cov_mode=0,
    random_state=None,
    tolerance=0.05,
):
    """
    Clusters sequences and splits data into train/test sets by grouping entire clusters.

    Parameters:
        df (pd.DataFrame): DataFrame with an ID column and a sequence column.
        sequence_col (str): Name of the column containing sequences.
        id_col (str): Name of the unique identifier column.
        test_size (float): Desired fraction of data in the test set (default 0.2).
        min_seq_id (float): Minimum sequence identity for clustering.
        coverage (float): Minimum alignment coverage for clustering.
        cov_mode (int): Coverage mode for clustering.
        random_state (int): Optional random state for reproducibility.
        tolerance (float): Acceptable deviation from test_size (default 0.05).

    Returns:
        (pd.DataFrame, pd.DataFrame): (train_df, test_df)
    """
    logger.info("Performing combined clustering and train/test split")
    logger.info(
        f"Parameters: sequence_col='{sequence_col}', id_col='{id_col}', test_size={test_size}"
    )
    logger.info(
        f"Clustering parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}"
    )

    _check_mmseqs()

    logger.info("Step 1: Clustering sequences")
    df_clustered = cluster(
        df=df,
        sequence_col=sequence_col,
        id_col=id_col,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode,
    )

    logger.info("Step 2: Splitting data based on sequence clusters")
    return split(
        df=df_clustered,
        group_col="representative_sequence",
        test_size=test_size,
        random_state=random_state,
        tolerance=tolerance,
    )


def train_test_val_cluster_split(
    df,
    sequence_col,
    id_col=None,
    test_size=0.2,
    val_size=0.1,
    min_seq_id=0.3,
    coverage=0.5,
    cov_mode=0,
    random_state=None,
    tolerance=0.05,
):
    """
    Clusters sequences and splits data into train, val, and test sets by grouping entire clusters.

    Parameters:
        df (pd.DataFrame): DataFrame with an ID column and a sequence column.
        sequence_col (str): Name of the column containing sequences.
        id_col (str): Name of the unique identifier column.
        test_size (float): Desired fraction of data in the test set (default 0.2).
        val_size (float): Desired fraction of data in the val set (default 0.1).
        min_seq_id (float): Minimum sequence identity for clustering.
        coverage (float): Minimum alignment coverage for clustering.
        cov_mode (int): Coverage mode for clustering.
        random_state (int): Optional random state for reproducibility.
        tolerance (float): Acceptable deviation from test_size and val_size (default 0.05).

    Returns:
        (pd.DataFrame, pd.DataFrame, pd.DataFrame): (train_df, val_df, test_df)
    """
    logger.info("Performing 3-way train/validation/test split with clustering")
    logger.info(f"Parameters: sequence_col='{sequence_col}', id_col='{id_col}'")
    logger.info(f"Split sizes: test_size={test_size}, val_size={val_size}")
    logger.info(
        f"Clustering parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}"
    )

    _check_mmseqs()

    logger.info("Step 1: Clustering sequences")
    df_clustered = cluster(
        df=df,
        sequence_col=sequence_col,
        id_col=id_col,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode,
    )

    logger.info("Step 2: Splitting into train+val vs test")
    train_val_df, test_df = split(
        df=df_clustered,
        group_col="representative_sequence",
        test_size=test_size,
        random_state=random_state,
        tolerance=tolerance,
    )

    logger.info("Step 3: Further splitting train+val into train vs val")
    adjusted_val_fraction = val_size / (1.0 - test_size)
    logger.debug(
        f"Adjusted validation fraction: {adjusted_val_fraction:.4f} of train+val set"
    )

    train_df, val_df = split(
        df=train_val_df,
        group_col="representative_sequence",
        test_size=adjusted_val_fraction,
        random_state=random_state,
        tolerance=tolerance,
    )

    total = len(df)
    logger.info("Final split results:")
    logger.info(f"  Train: {len(train_df)} sequences ({len(train_df)/total:.2%})")
    logger.info(f"  Validation: {len(val_df)} sequences ({len(val_df)/total:.2%})")
    logger.info(f"  Test: {len(test_df)} sequences ({len(test_df)/total:.2%})")

    return train_df, val_df, test_df


def constrained_split(
    df,
    group_col="representative_sequence",
    id_col=None,
    test_size=0.2,
    random_state=None,
    tolerance=0.05,
    force_train_ids=None,
    force_test_ids=None,
    id_type="sequence",
):
    """
    Splits data with constraints on which sequences or groups must be in the train or test set.
    Assumes that the DataFrame already has a grouping column (e.g., 'representative_sequence') that defines clusters.

    Parameters:
        df (pd.DataFrame): DataFrame containing sequences.
        group_col (str): Column indicating the group/cluster each sequence belongs to.
        id_col (str): Name of the unique identifier column (used when id_type is "sequence").
        test_size (float): Desired fraction of data in the test set (default 0.2).
        random_state (int): Random seed for reproducibility.
        tolerance (float): Acceptable deviation from test_size (default 0.05).
        force_train_ids (list): IDs (or groups, depending on id_type) that must be in the training set.
        force_test_ids (list): IDs (or groups, depending on id_type) that must be in the test set.
        id_type (str): Specifies whether the forced IDs refer to individual sequences ("sequence", uses id_col)
                       or groups ("cluster" or "group", uses group_col). Default is "sequence".

    Returns:
        (pd.DataFrame, pd.DataFrame): Tuple of (train_df, test_df)
    """
    logger.info("Performing constrained train/test split")

    # Set the random seed if provided
    if random_state is not None:
        set_seed(random_state)

    # Initialize forced IDs if not provided
    force_train_ids = [] if force_train_ids is None else force_train_ids
    force_test_ids = [] if force_test_ids is None else force_test_ids

    # Check for conflicts at the sequence level
    conflicts = set(force_train_ids).intersection(set(force_test_ids))
    if conflicts:
        logger.error(
            f"Found {len(conflicts)} IDs in both force_train_ids and force_test_ids"
        )
        raise ValueError("Cannot force the same IDs to both train and test sets")

    forced_train_groups = set()
    forced_test_groups = set()

    if id_type == "sequence":
        if id_col is None:
            raise ValueError("id_col must be provided when id_type is 'sequence'")
        # Determine forced groups based on forced sequence IDs
        if force_train_ids:
            train_mask = df[id_col].isin(force_train_ids)
            forced_train_groups.update(df.loc[train_mask, group_col].unique())
            logger.info(
                f"Forcing {len(forced_train_groups)} groups to train based on {len(force_train_ids)} sequence IDs"
            )
        if force_test_ids:
            test_mask = df[id_col].isin(force_test_ids)
            forced_test_groups.update(df.loc[test_mask, group_col].unique())
            logger.info(
                f"Forcing {len(forced_test_groups)} groups to test based on {len(force_test_ids)} sequence IDs"
            )
    elif id_type in ["cluster", "group"]:
        forced_train_groups.update(force_train_ids)
        forced_test_groups.update(force_test_ids)
        logger.info(
            f"Forcing {len(forced_train_groups)} groups to train and {len(forced_test_groups)} to test"
        )
    else:
        raise ValueError(f"Invalid id_type: {id_type}. Must be 'sequence' or 'cluster'")

    # Check for conflicts at the group level
    group_conflicts = forced_train_groups.intersection(forced_test_groups)
    if group_conflicts:
        logger.error(
            f"Found {len(group_conflicts)} groups forced to both train and test"
        )
        raise ValueError(
            "Constraint conflict: some groups are forced to both train and test sets"
        )

    # Create forced splits based on group_col membership
    train_forced_mask = df[group_col].isin(forced_train_groups)
    test_forced_mask = df[group_col].isin(forced_test_groups)

    train_forced = df[train_forced_mask]
    test_forced = df[test_forced_mask]
    remaining = df[~(train_forced_mask | test_forced_mask)]

    total_size = len(df)
    train_forced_size = len(train_forced)
    test_forced_size = len(test_forced)
    remaining_size = len(remaining)

    logger.info(
        f"Pre-assigned {train_forced_size} sequences to train ({train_forced_size/total_size:.2%})"
    )
    logger.info(
        f"Pre-assigned {test_forced_size} sequences to test ({test_forced_size/total_size:.2%})"
    )
    logger.info(
        f"Remaining {remaining_size} sequences to split ({remaining_size/total_size:.2%})"
    )

    if test_forced_size / total_size > test_size + tolerance:
        logger.warning(
            f"Forced test assignments ({test_forced_size/total_size:.2%}) exceed desired test size ({test_size:.2%}) by more than tolerance ({tolerance:.2%})"
        )

    # Calculate target test size for the remaining sequences
    target_test_from_remaining = max(0, (test_size * total_size) - test_forced_size)
    if remaining_size > 0:
        adjusted_test_size = target_test_from_remaining / remaining_size
        adjusted_test_size = min(1.0, max(0.0, adjusted_test_size))
        logger.info(f"Adjusted test size for remaining data: {adjusted_test_size:.2%}")

        # Split remaining data using the existing group-aware split function
        if 0 < adjusted_test_size < 1:
            train_remaining, test_remaining = split(
                df=remaining,
                group_col=group_col,
                test_size=adjusted_test_size,
                random_state=random_state,
                tolerance=tolerance,
            )
        elif adjusted_test_size <= 0:
            train_remaining, test_remaining = remaining.copy(), remaining.iloc[0:0]
            logger.info("All remaining sequences assigned to train set")
        else:  # adjusted_test_size >= 1
            train_remaining, test_remaining = remaining.iloc[0:0], remaining.copy()
            logger.info("All remaining sequences assigned to test set")
    else:
        train_remaining = pd.DataFrame(columns=df.columns)
        test_remaining = pd.DataFrame(columns=df.columns)

    # Combine forced and split remaining data
    train_df = pd.concat([train_forced, train_remaining])
    test_df = pd.concat([test_forced, test_remaining])

    final_train_pct = len(train_df) / total_size
    final_test_pct = len(test_df) / total_size

    logger.info("Final split results:")
    logger.info(f"  Train: {len(train_df)} sequences ({final_train_pct:.2%})")
    logger.info(f"  Test: {len(test_df)} sequences ({final_test_pct:.2%})")

    if abs(final_test_pct - test_size) > tolerance:
        logger.warning(
            f"Final test fraction ({final_test_pct:.2%}) differs from target ({test_size:.2%}) by more than tolerance ({tolerance:.2%})"
        )

    return train_df, test_df


def cluster_kfold(
    df,
    sequence_col,
    id_col=None,
    n_splits=5,
    min_seq_id=0.3,
    coverage=0.5,
    cov_mode=0,
    alignment_mode=0,
    random_state=None,
    shuffle=True,
    return_indices=False,
):
    """
    Performs k-fold cross-validation while respecting sequence clustering.
    Ensures that sequences in the same cluster are always in the same fold.

    Parameters:
        df (pd.DataFrame): DataFrame with sequences
        sequence_col (str): Name of the column containing sequences
        id_col (str): Name of the unique identifier column
        n_splits (int): Number of folds for cross-validation
        min_seq_id (float): Minimum sequence identity for clustering
        coverage (float): Minimum alignment coverage for clustering
        cov_mode (int): Coverage mode for clustering
        alignment_mode (int): Alignment mode for clustering
        random_state (int): Random seed for reproducibility
        shuffle (bool): Whether to shuffle clusters before assigning to folds
        return_indices (bool): If True, returns indices instead of DataFrames

    Returns:
        list: List of (train_df, test_df) tuples for each fold
              If return_indices=True, returns (train_indices, test_indices) instead
    """
    logger.info(
        f"Performing {n_splits}-fold cross-validation with cluster-aware splits"
    )

    # Set the random seed if provided
    if random_state is not None:
        set_seed(random_state)

    # First, cluster the sequences
    df_clustered = cluster(
        df=df,
        sequence_col=sequence_col,
        id_col=id_col,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode,
        alignment_mode=alignment_mode,
    )

    # Get unique clusters and their sizes
    cluster_sizes = df_clustered.groupby("representative_sequence").size()
    clusters = list(cluster_sizes.index)

    logger.info(f"Found {len(clusters)} unique sequence clusters")

    if len(clusters) < n_splits:
        logger.warning(
            f"Number of clusters ({len(clusters)}) is less than requested folds ({n_splits}). "
            f"Some folds will be empty."
        )

    # Sort clusters by size (descending) for better balancing
    clusters_with_sizes = [(cluster, cluster_sizes[cluster]) for cluster in clusters]
    clusters_with_sizes.sort(key=lambda x: x[1], reverse=True)

    # Optional shuffling of similarly-sized clusters to ensure randomness
    if shuffle:
        # Group clusters by size
        size_groups = {}
        for cluster, size in clusters_with_sizes:
            if size not in size_groups:
                size_groups[size] = []
            size_groups[size].append(cluster)

        # Shuffle each size group
        for size in size_groups:
            random.shuffle(size_groups[size])

        # Reconstruct clusters_with_sizes with preserved size ordering but shuffled within size
        clusters_with_sizes = []
        for size in sorted(size_groups.keys(), reverse=True):
            for cluster in size_groups[size]:
                clusters_with_sizes.append((cluster, size))

    # Initialize fold assignments
    fold_sizes = [0] * n_splits
    cluster_to_fold = {}

    # Assign each cluster to the smallest fold
    for cluster, size in clusters_with_sizes:
        smallest_fold = fold_sizes.index(min(fold_sizes))
        fold_sizes[smallest_fold] += size
        cluster_to_fold[cluster] = smallest_fold

    # Report fold balance
    total_size = sum(fold_sizes)
    fold_percentages = [size / total_size for size in fold_sizes]

    logger.info(f"Fold sizes: {fold_sizes}")
    logger.info(f"Fold percentages: {[f'{p:.1%}' for p in fold_percentages]}")

    # Check for imbalanced folds
    fold_imbalance = max(fold_percentages) - min(fold_percentages)
    if fold_imbalance > 0.1:  # More than 10% difference
        logger.warning(
            f"Folds are imbalanced (max={max(fold_percentages):.1%}, min={min(fold_percentages):.1%}, "
            f"diff={fold_imbalance:.1%}). Consider adjusting clustering parameters."
        )

    # For each fold, create train/test split
    result = []
    for fold_idx in range(n_splits):
        logger.info(f"Preparing fold {fold_idx+1}/{n_splits}")

        # Create mask for test set (current fold)
        test_mask = df_clustered["representative_sequence"].apply(
            lambda x: cluster_to_fold.get(x) == fold_idx
        )

        if return_indices:
            train_indices = df_clustered[~test_mask].index
            test_indices = df_clustered[test_mask].index
            result.append((train_indices, test_indices))
        else:
            train_df = df_clustered[~test_mask].copy()
            test_df = df_clustered[test_mask].copy()
            result.append((train_df, test_df))

        n_train = len(train_indices if return_indices else train_df)
        n_test = len(test_indices if return_indices else test_df)
        logger.info(
            f"Fold {fold_idx+1}: train={n_train} ({n_train/(n_train+n_test):.1%}), "
            f"test={n_test} ({n_test/(n_train+n_test):.1%})"
        )

    return result
