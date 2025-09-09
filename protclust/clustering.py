import logging
import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .logger import logger
from .utils import _check_mmseqs, check_random_state


def cluster(
    df: pd.DataFrame,
    sequence_col: str,
    id_col: Optional[str] = None,
    min_seq_id: float = 0.3,
    coverage: float = 0.8,
    cov_mode: int = 0,
    alignment_mode: int = 0,
    cluster_mode: int = 0,
    cluster_steps: int = 1,
    threads: int = 1,
    random_state: Optional[Union[int, np.random.RandomState]] = None,
    temp_dir: Optional[str] = None,
    return_stats: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, Any]]]:
    """
    Clusters protein sequences using MMseqs2.

    Parameters:
        df (pd.DataFrame): Input DataFrame with sequences.
        sequence_col (str): Name of the column containing sequences.
        id_col (str, optional): ID column name. Uses index if None.
        min_seq_id (float): Minimum sequence identity (0.0-1.0, default 0.3).
        coverage (float): Minimum coverage (0.0-1.0, default 0.8).
        cov_mode (int): Coverage mode (0-2, default 0):
            0: coverage of query and target
            1: coverage of target
            2: coverage of query
        alignment_mode (int): Alignment mode (0-4, default 0):
            0: automatic
            1: only score
            2: only ungapped alignment
            3: score and end_pos
            4: ungapped alignment and end_pos
        cluster_mode (int): Clustering mode (0-2, default 0):
            0: Set-Cover (greedy)
            1: Connected component (BLASTclust)
            2: Greedy by sequence length (CD-HIT)
        cluster_steps (int): Number of cascaded clustering steps (default 1).
        threads (int): Number of threads (default 1).
        random_state (int or RandomState, optional): For reproducibility.
        temp_dir (str, optional): Directory for temporary files.
        return_stats (bool): Return clustering statistics.

    Returns:
        DataFrame with 'cluster_representative' column, optionally with stats.
    """
    logger.info("Starting sequence clustering with MMseqs2")
    logger.info(
        f"Parameters: min_seq_id={min_seq_id}, coverage={coverage}, cov_mode={cov_mode}, alignment_mode={alignment_mode}, cluster_mode={cluster_mode}, cluster_steps={cluster_steps}"
    )

    # Validation
    _check_mmseqs()

    if not 0.0 <= min_seq_id <= 1.0:
        raise ValueError(f"min_seq_id must be between 0.0 and 1.0, got {min_seq_id}")
    if not 0.0 <= coverage <= 1.0:
        raise ValueError(f"coverage must be between 0.0 and 1.0, got {coverage}")
    if cov_mode not in [0, 1, 2]:
        raise ValueError(f"cov_mode must be 0, 1, or 2, got {cov_mode}")
    if alignment_mode not in [0, 1, 2, 3, 4]:
        raise ValueError(f"alignment_mode must be 0-4, got {alignment_mode}")
    if cluster_mode not in [0, 1, 2]:
        raise ValueError(f"cluster_mode must be 0, 1, or 2, got {cluster_mode}")
    if not isinstance(cluster_steps, int) or cluster_steps <= 0:
        raise ValueError(f"cluster_steps must be a positive integer, got {cluster_steps}")
    if threads < 1:
        raise ValueError(f"threads must be >= 1, got {threads}")

    # Get random state
    random_state = check_random_state(random_state)

    # Create a deep copy to avoid SettingWithCopyWarning
    result_df = df.copy(deep=True)

    if id_col is None:
        result_df = result_df.reset_index()
        id_col = "index"
        logger.debug(f"No id_col provided, using '{id_col}' as identifier")

    if sequence_col not in result_df or id_col not in result_df:
        logger.error(f"Required columns missing: {sequence_col} or {id_col}")
        raise ValueError(f"The DataFrame must have '{id_col}' and '{sequence_col}'.")

    logger.info(f"Clustering {len(result_df)} sequences")

    # Use .loc for assignment to avoid SettingWithCopyWarning
    result_df.loc[:, "sanitized_id"] = result_df[id_col].str.replace(" ", "_")

    # If user specified random state, sort input sequences to ensure reproducibility
    if random_state is not None:
        result_df = result_df.sort_values(by=[sequence_col, "sanitized_id"])

    # Create or use the specified temporary directory
    if temp_dir is None:
        tmp_dir = tempfile.mkdtemp()
        created_temp_dir = True
    else:
        os.makedirs(temp_dir, exist_ok=True)
        tmp_dir = temp_dir
        created_temp_dir = False

    logger.debug(f"Using temporary directory: {tmp_dir}")

    # Dictionary to store clustering statistics
    clustering_stats = {
        "parameters": {
            "min_seq_id": min_seq_id,
            "coverage": coverage,
            "cov_mode": cov_mode,
            "alignment_mode": alignment_mode,
            "threads": threads,
            "random_state": random_state.randint(0, 2**32 - 1)
            if random_state is not None
            else None,
        },
        "input_sequences": len(result_df),
    }

    try:
        input_fasta = os.path.join(tmp_dir, "input.fasta")
        with open(input_fasta, "w") as fasta_file:
            for _, row in result_df.iterrows():
                fasta_file.write(f">{row['sanitized_id']}\n{row[sequence_col]}\n")

        logger.debug(f"Wrote {len(result_df)} sequences to FASTA file")

        # Always use the full pipeline for consistency
        clusters_file = _run_clustering_pipeline(
            input_fasta=input_fasta,
            tmp_dir=tmp_dir,
            min_seq_id=min_seq_id,
            coverage=coverage,
            cov_mode=cov_mode,
            alignment_mode=alignment_mode,
            cluster_mode=cluster_mode,
            cluster_steps=cluster_steps,
            threads=threads,
            random_state=random_state,
        )

        if not os.path.exists(clusters_file):
            logger.error("MMseqs2 clustering results file not found")
            raise FileNotFoundError("MMseqs2 clustering results not found.")

        logger.debug(f"Reading clustering results from {clusters_file}")

        # Parse clustering results
        cluster_map = {}
        cluster_sizes = {}
        with open(clusters_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    rep, seq = parts[0], parts[1]
                    cluster_map[seq] = rep
                    cluster_sizes[rep] = cluster_sizes.get(rep, 0) + 1

        logger.info(f"Found {len(cluster_sizes)} clusters")

        # Calculate cluster statistics
        clustering_stats["num_clusters"] = len(cluster_sizes)
        clustering_stats["singleton_clusters"] = sum(
            1 for size in cluster_sizes.values() if size == 1
        )
        clustering_stats["largest_cluster_size"] = (
            max(cluster_sizes.values()) if cluster_sizes else 0
        )
        clustering_stats["average_cluster_size"] = (
            np.mean(list(cluster_sizes.values())) if cluster_sizes else 0
        )

        # Detailed cluster size distribution
        cluster_size_counts = {}
        for size in cluster_sizes.values():
            cluster_size_counts[size] = cluster_size_counts.get(size, 0) + 1

        clustering_stats["cluster_size_distribution"] = cluster_size_counts

        # Create mapping from sanitized IDs back to original IDs
        reverse_map = dict(zip(result_df["sanitized_id"], result_df[id_col]))

        # Use .loc for assignment to avoid SettingWithCopyWarning
        result_df.loc[:, "cluster_representative"] = result_df["sanitized_id"].apply(
            lambda x: reverse_map.get(cluster_map.get(x, x), x)
        )

        logger.info("Clustering complete, added 'cluster_representative' column to DataFrame")

    finally:
        # Clean up temporary directory
        if created_temp_dir:
            logger.debug(f"Cleaning up temporary directory: {tmp_dir}")
            shutil.rmtree(tmp_dir)

    # Avoid inplace=True to prevent SettingWithCopyWarning
    result_df = result_df.drop(columns=["sanitized_id"])

    # Return results based on return_stats flag
    if return_stats:
        return result_df, clustering_stats
    else:
        return result_df


def _run_clustering_pipeline(
    input_fasta: str,
    tmp_dir: str,
    min_seq_id: float,
    coverage: float,
    cov_mode: int,
    alignment_mode: int,
    cluster_mode: int,
    cluster_steps: int,
    threads: int,
    random_state: Optional[np.random.RandomState] = None,
) -> str:
    """Run MMseqs2 clustering pipeline, choosing optimal method based on sequence identity."""

    # Use linclust for high identity (faster), cluster for lower identity
    use_linclust = min_seq_id >= 0.8

    # Define paths
    db_path = os.path.join(tmp_dir, "seqDB")
    cluster_db_path = os.path.join(tmp_dir, "clusterDB")
    cluster_tsv_path = os.path.join(tmp_dir, "cluster.tsv")
    tmp_mmseqs = os.path.join(tmp_dir, "tmp_mmseqs")

    os.makedirs(tmp_mmseqs, exist_ok=True)

    # Step 1: Create sequence database
    createdb_cmd = ["mmseqs", "createdb", input_fasta, db_path]
    logger.debug(f"Running MMseqs2 createdb command: {' '.join(createdb_cmd)}")

    try:
        subprocess.run(createdb_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"MMseqs2 createdb failed: {e.stderr}")
        raise RuntimeError(
            "Failed to create MMseqs2 database. Check that input sequences are valid."
        ) from e

    # Step 2: Run clustering
    if use_linclust:
        # Use linclust for high sequence identity
        cluster_cmd = [
            "mmseqs",
            "linclust",
            db_path,
            cluster_db_path,
            tmp_mmseqs,
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
        ]
    else:
        # Use cluster for lower sequence identity
        cluster_cmd = [
            "mmseqs",
            "cluster",
            db_path,
            cluster_db_path,
            tmp_mmseqs,
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
        ]

        # Only add alignment_mode for cluster (not supported by linclust)
        if alignment_mode != 0:
            cluster_cmd.extend(["--alignment-mode", str(alignment_mode)])

    # Add cluster mode (supported by both linclust and cluster)
    if cluster_mode != 0:
        cluster_cmd.extend(["--cluster-mode", str(cluster_mode)])

    # Add cluster steps for cascaded clustering
    if cluster_steps > 1:
        cluster_cmd.extend(["--cluster-steps", str(cluster_steps)])

    # Handle threading and reproducibility
    if random_state is not None:
        cluster_cmd.extend(["--threads", "1"])
        if threads > 1:
            logger.warning("Setting threads=1 to ensure reproducibility with random_state")
    elif threads > 1:
        cluster_cmd.extend(["--threads", str(threads)])

    logger.debug(
        f"Running MMseqs2 {'linclust' if use_linclust else 'cluster'} command: {' '.join(cluster_cmd)}"
    )

    try:
        if logger.level <= logging.DEBUG:
            subprocess.run(cluster_cmd, check=True, text=True)
        else:
            subprocess.run(
                cluster_cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        logger.error(f"MMseqs2 clustering failed: {e.stderr}")
        raise RuntimeError(
            "MMseqs2 clustering failed. Check clustering parameters and sequence format."
        ) from e

    # Step 3: Create cluster TSV file
    createtsv_cmd = ["mmseqs", "createtsv", db_path, db_path, cluster_db_path, cluster_tsv_path]

    logger.debug(f"Running MMseqs2 createtsv command: {' '.join(createtsv_cmd)}")

    try:
        subprocess.run(createtsv_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"MMseqs2 createtsv failed: {e.stderr}")
        raise RuntimeError("Failed to create cluster results file.") from e

    return cluster_tsv_path
