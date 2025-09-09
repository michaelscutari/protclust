"""Tests for split integrity using synthetic data with predictable patterns."""

import os
import subprocess
import tempfile

import pandas as pd

from protclust import cluster, split
from protclust.utils import _check_mmseqs


def test_basic_split_integrity(synthetic_cluster_data, mmseqs_installed):
    """Test that splits maintain cluster integrity."""
    df = synthetic_cluster_data.copy()

    # First, cluster the data
    clustered_df = cluster(
        df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.8,
        coverage=0.8,
    )

    # Perform splitting with different test sizes
    test_sizes = [0.2, 0.3, 0.4]

    for test_size in test_sizes:
        # Split the data
        train_df, test_df = split(
            clustered_df,
            group_col="cluster_representative",
            test_size=test_size,
            random_state=42,
        )

        # Check that all rows are accounted for
        assert len(train_df) + len(test_df) == len(clustered_df)

        # Check that clusters are kept intact (no cluster spans both train and test)
        train_clusters = set(train_df["cluster_representative"])
        test_clusters = set(test_df["cluster_representative"])

        assert len(train_clusters.intersection(test_clusters)) == 0, (
            f"Split with test_size={test_size} has clusters that span both train and test sets"
        )

        # Check that the test size is reasonably close to target
        actual_test_size = len(test_df) / len(clustered_df)
        assert abs(actual_test_size - test_size) <= 0.1, (
            f"Target test_size={test_size}, but got {actual_test_size:.2f}"
        )


def test_split_reproducibility(synthetic_cluster_data, mmseqs_installed):
    """Test that splits are reproducible with the same random seed."""
    df = synthetic_cluster_data.copy()

    # First, cluster the data
    clustered_df = cluster(
        df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.8,
        coverage=0.8,
    )

    # Perform the split twice with the same seed
    train1, test1 = split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        random_state=42,
    )

    train2, test2 = split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        random_state=42,
    )

    # Check that the splits are identical
    assert set(train1["id"]) == set(train2["id"]), "Train sets differ despite same random seed"
    assert set(test1["id"]) == set(test2["id"]), "Test sets differ despite same random seed"

    # Perform another split with a different seed
    train3, test3 = split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        random_state=24,
    )

    # Check that the splits are different
    assert set(train1["id"]) != set(train3["id"]), (
        "Train sets are identical despite different seeds"
    )
    assert set(test1["id"]) != set(test3["id"]), "Test sets are identical despite different seeds"


def test_split_with_property_distribution(synthetic_cluster_data, mmseqs_installed):
    """Test that property distributions are maintained in splits."""
    df = synthetic_cluster_data.copy()

    # First, cluster the data
    clustered_df = cluster(
        df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.8,
        coverage=0.8,
    )

    # Split the data
    train_df, test_df = split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        random_state=42,
    )

    # Check property value distributions
    # Since our synthetic data correlates property_value with cluster_id,
    # we can check if the property distribution is roughly maintained

    # Get overall stats
    overall_mean = df["property_value"].mean()
    overall_std = df["property_value"].std()

    # Get split stats
    train_mean = train_df["property_value"].mean()
    test_mean = test_df["property_value"].mean()

    # Check that means are reasonably close
    # Allow for some deviation due to cluster-based splitting
    assert abs(train_mean - overall_mean) < overall_std, (
        f"Train set property mean ({train_mean:.2f}) differs too much from overall ({overall_mean:.2f})"
    )
    assert abs(test_mean - overall_mean) < overall_std, (
        f"Test set property mean ({test_mean:.2f}) differs too much from overall ({overall_mean:.2f})"
    )


def test_sequence_identity_between_splits_mmseqs(synthetic_cluster_data, mmseqs_installed):
    """
    Test that sequence identity between train and test sets remains below the specified threshold.
    This test uses MMseqs2 directly to verify sequence identity, ensuring the same algorithm
    that's used for clustering is also used for verification.
    """
    _check_mmseqs()  # Ensure MMseqs2 is installed
    df = synthetic_cluster_data.copy()

    # Test different clustering thresholds
    min_seq_id_values = [0.5, 0.7, 0.9]

    for min_seq_id in min_seq_id_values:
        # First cluster, then split
        clustered_df = cluster(
            df,
            sequence_col="sequence",
            id_col="id",
            min_seq_id=min_seq_id,  # Vary the sequence identity threshold
            coverage=0.8,
            random_state=42,
        )

        train_df, test_df = split(
            clustered_df,
            group_col="cluster_representative",
            test_size=0.3,
            random_state=42,
        )

        # Skip if either split is too small
        if len(train_df) < 2 or len(test_df) < 2:
            continue

        # Create temporary files and directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            train_fasta = os.path.join(tmp_dir, "train.fasta")
            test_fasta = os.path.join(tmp_dir, "test.fasta")
            results_file = os.path.join(tmp_dir, "results.tsv")

            # Write train sequences to FASTA
            with open(train_fasta, "w") as f:
                for _, row in train_df.iterrows():
                    f.write(f">train_{row['id']}\n{row['sequence']}\n")

            # Write test sequences to FASTA
            with open(test_fasta, "w") as f:
                for _, row in test_df.iterrows():
                    f.write(f">test_{row['id']}\n{row['sequence']}\n")

            # Create MMseqs2 databases
            train_db = os.path.join(tmp_dir, "trainDB")
            test_db = os.path.join(tmp_dir, "testDB")
            result_db = os.path.join(tmp_dir, "resultDB")
            tmp_mmseqs = os.path.join(tmp_dir, "tmp")

            # Convert FASTA to MMseqs2 DB format
            subprocess.run(["mmseqs", "createdb", train_fasta, train_db], check=True)
            subprocess.run(["mmseqs", "createdb", test_fasta, test_db], check=True)

            # Run search to find sequence similarities
            # Set e-value threshold high to capture all possible matches
            subprocess.run(
                [
                    "mmseqs",
                    "search",
                    train_db,
                    test_db,
                    result_db,
                    tmp_mmseqs,
                    "--threads",
                    "1",  # For reproducibility
                    "-e",
                    "100",  # High e-value to find all matches
                    "--cov-mode",
                    "0",  # Coverage of query and target
                    "-c",
                    "0.0",  # Any coverage is fine for this test
                ],
                check=True,
            )

            # Convert results to TSV format with sequence identity
            subprocess.run(
                [
                    "mmseqs",
                    "convertalis",
                    train_db,
                    test_db,
                    result_db,
                    results_file,
                    "--format-output",
                    "query,target,pident,fident",
                ],
                check=True,
            )

            # Load and analyze results
            if os.path.exists(results_file) and os.path.getsize(results_file) > 0:
                results_df = pd.read_csv(results_file, sep="\t", header=None)
                results_df.columns = ["query", "target", "pident", "fident"]

                # pident is percentage identity (0-100), convert to fraction (0-1)
                max_identity = results_df["pident"].max() / 100.0 if not results_df.empty else 0.0

                # Check that maximum identity is below threshold (with small tolerance)
                tolerance = 0.05  # Allow 5% tolerance
                assert max_identity <= min_seq_id + tolerance, (
                    f"With min_seq_id={min_seq_id}, found sequences between splits "
                    f"with identity {max_identity:.2f}, which exceeds the threshold"
                )
