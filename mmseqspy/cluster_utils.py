import os
import tempfile
import subprocess
import shutil
import pandas as pd

def _check_mmseqs():
    """
    Ensures 'mmseqs' command is in PATH.
    """
    if shutil.which("mmseqs") is None:
        raise EnvironmentError(
            "MMseqs2 is not installed or not found in PATH. "
            "See the README for installation instructions."
        )

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
    df[sequence_col] = df[sequence_col].str.upper()
    df = df.dropna(subset=[sequence_col])
    valid_sequence_mask = df[sequence_col].apply(
        lambda seq: all(aa in valid_amino_acids for aa in seq)
    )
    return df[valid_sequence_mask].reset_index(drop=True)

def cluster(
    df,
    sequence_col,
    id_col="id",
    min_seq_id=0.3,
    coverage=0.5,
    cov_mode=0
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

    Returns:
        pd.DataFrame: Original DataFrame with a new 'representative_sequence' column.
    """
    _check_mmseqs()

    if sequence_col not in df or id_col not in df:
        raise ValueError(f"The DataFrame must have '{id_col}' and '{sequence_col}'.")

    df["sanitized_id"] = df[id_col].str.replace(" ", "_")

    tmp_dir = tempfile.mkdtemp()
    try:
        input_fasta = os.path.join(tmp_dir, "input.fasta")
        with open(input_fasta, "w") as fasta_file:
            for _, row in df.iterrows():
                fasta_file.write(f">{row['sanitized_id']}\n{row[sequence_col]}\n")

        output_dir = os.path.join(tmp_dir, "output")
        tmp_mmseqs = os.path.join(tmp_dir, "tmp_mmseqs")
        subprocess.run([
            "mmseqs", "easy-cluster", input_fasta, output_dir, tmp_mmseqs,
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode)
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        clusters_file = os.path.join(output_dir + "_cluster.tsv")
        if not os.path.exists(clusters_file):
            raise FileNotFoundError("MMseqs2 clustering results not found.")

        cluster_map = {}
        with open(clusters_file, "r") as f:
            for line in f:
                rep, seq = line.strip().split("\t")
                cluster_map[seq] = rep

        reverse_map = dict(zip(df["sanitized_id"], df[id_col]))
        df["representative_sequence"] = df["sanitized_id"].apply(
            lambda x: reverse_map.get(cluster_map.get(x, x), x)
        )

    finally:
        shutil.rmtree(tmp_dir)

    df.drop(columns=["sanitized_id"], inplace=True)
    return df

def split(
    df,
    group_col,
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
    # Determine the total number of rows
    total_sequences = len(df)
    target_test_count = int(round(test_size * total_sequences))

    # Compute group sizes
    group_sizes_df = df.groupby(group_col).size().reset_index(name="group_size")
    groups = group_sizes_df[group_col].tolist()
    sizes = group_sizes_df["group_size"].tolist()

    # Subset-sum dynamic programming to find the best combination of groups
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

    # Split into train and test based on chosen groups
    test_df = df[df[group_col].isin(chosen_groups)]
    train_df = df[~df[group_col].isin(chosen_groups)]

    # Check how close we are to the desired fraction
    achieved_test_fraction = len(test_df) / total_sequences
    if abs(achieved_test_fraction - test_size) > tolerance:
        print(
            f"Warning: Desired test fraction = {test_size:.2f}, "
            f"achieved = {achieved_test_fraction:.2f}. "
            "Closest possible split."
        )

    return train_df, test_df

def cluster_split(
    df,
    sequence_col,
    id_col="id",
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
    _check_mmseqs()
    # Cluster the data to identify representative sequences
    df_clustered = cluster(
        df=df,
        sequence_col=sequence_col,
        id_col=id_col,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode
    )
    # Perform the split by grouping on 'representative_sequence'
    return split(
        df=df_clustered,
        group_col="representative_sequence",
        test_size=test_size,
        random_state=random_state,
        tolerance=tolerance
    )