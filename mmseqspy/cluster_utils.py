import os
import tempfile
import subprocess
import shutil
import pandas as pd
import logging

# Set up logging
logger = logging.getLogger("protein_clustering")
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)  # Default level

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

def clean(
    df, 
    sequence_col='sequence', 
    valid_amino_acids='ACDEFGHIKLMNPQRSTVWY'
):
    """
    Removes sequences with invalid protein characters.

    Parameters:
        df (pd.DataFrame): Input DataFrame with protein sequences.
        sequence_col (str): Name of the column containing sequences.
        valid_amino_acids (str): String of valid amino acid characters.

    Returns:
        pd.DataFrame: Cleaned DataFrame with only valid sequences.
    """
    logger.info(f"Cleaning sequences in column '{sequence_col}' with valid amino acids: {valid_amino_acids}")
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
    logger.info(f"Parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}, alignment_mode={alignment_mode}")
    
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
            "mmseqs", "easy-cluster", input_fasta, output_dir, tmp_mmseqs,
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
            '--alignment-mode', str(alignment_mode)
        ]
        
        logger.debug(f"Running MMseqs2 command: {' '.join(mmseqs_cmd)}")
        
        if logger.level <= logging.DEBUG:
            subprocess.run(mmseqs_cmd, check=True)
        else:
            subprocess.run(mmseqs_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
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
        
        logger.info("Clustering complete, added 'representative_sequence' column to DataFrame")
        
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
    tolerance=0.05
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
    
    logger.info(f"Train set: {len(train_df)} sequences ({len(train_df) / total_sequences:.2%})")
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
    tolerance=0.05
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
    logger.info(f"Parameters: sequence_col='{sequence_col}', id_col='{id_col}', test_size={test_size}")
    logger.info(f"Clustering parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}")
    
    _check_mmseqs()
    
    logger.info("Step 1: Clustering sequences")
    df_clustered = cluster(
        df=df,
        sequence_col=sequence_col,
        id_col=id_col,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode
    )
    
    logger.info("Step 2: Splitting data based on sequence clusters")
    return split(
        df=df_clustered,
        group_col="representative_sequence",
        test_size=test_size,
        random_state=random_state,
        tolerance=tolerance
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
    tolerance=0.05
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
    logger.info(f"Clustering parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}")
    
    _check_mmseqs()
    
    logger.info("Step 1: Clustering sequences")
    df_clustered = cluster(
        df=df,
        sequence_col=sequence_col,
        id_col=id_col,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode
    )
    
    logger.info("Step 2: Splitting into train+val vs test")
    train_val_df, test_df = split(
        df=df_clustered,
        group_col="representative_sequence",
        test_size=test_size,
        random_state=random_state,
        tolerance=tolerance
    )
    
    logger.info("Step 3: Further splitting train+val into train vs val")
    adjusted_val_fraction = val_size / (1.0 - test_size)
    logger.debug(f"Adjusted validation fraction: {adjusted_val_fraction:.4f} of train+val set")
    
    train_df, val_df = split(
        df=train_val_df,
        group_col="representative_sequence",
        test_size=adjusted_val_fraction,
        random_state=random_state,
        tolerance=tolerance
    )
    
    total = len(df)
    logger.info(f"Final split results:")
    logger.info(f"  Train: {len(train_df)} sequences ({len(train_df)/total:.2%})")
    logger.info(f"  Validation: {len(val_df)} sequences ({len(val_df)/total:.2%})")
    logger.info(f"  Test: {len(test_df)} sequences ({len(test_df)/total:.2%})")
    
    return train_df, val_df, test_df