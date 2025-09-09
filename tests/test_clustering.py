from protclust import cluster


def test_cluster_sequences(fluorescence_data, mmseqs_installed):
    """Test clustering protein sequences."""
    # Make a copy of data
    df = fluorescence_data.copy()

    # Run clustering
    clustered_df = cluster(
        df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.5,  # 50% sequence identity threshold
        coverage=0.8,  # 80% coverage
    )

    # Check results
    assert "cluster_representative" in clustered_df.columns
    assert len(clustered_df) == len(df)  # Should preserve all rows

    # Count unique clusters
    n_clusters = clustered_df["cluster_representative"].nunique()

    # Basic sanity check - clusters should be fewer than sequences
    assert 1 <= n_clusters <= len(df)

    # Check that all representative_sequence values exist in the id column
    assert set(clustered_df["cluster_representative"]).issubset(set(clustered_df["id"]))


def test_cluster_debug_logging(fluorescence_data, mmseqs_installed, monkeypatch):
    """Test clustering with debug logging to cover verbose output paths."""
    import logging

    from protclust.logger import logger

    # Store original level and set up logging capture
    original_level = logger.level
    logger.setLevel(logging.DEBUG)

    try:
        # Run clustering with a small dataset
        df = fluorescence_data.head(10).copy()
        result = cluster(df, sequence_col="sequence", id_col="id")

        # Verify the clustering completed successfully
        assert "cluster_representative" in result.columns

    finally:
        # Restore original logging level
        logger.setLevel(original_level)


def test_cluster_with_debug_output(fluorescence_data, mmseqs_installed, caplog):
    """Test clustering with debug output enabled."""
    import logging

    from protclust import cluster
    from protclust.logger import logger

    # Set debug level temporarily
    original_level = logger.level
    logger.setLevel(logging.DEBUG)

    try:
        # Just a small sample to keep it fast
        df_small = fluorescence_data.head(10).copy()
        result = cluster(df_small, sequence_col="sequence", id_col="id")

        # Verify it ran successfully
        assert "cluster_representative" in result.columns
    finally:
        # Restore original level
        logger.setLevel(original_level)


def test_cluster_reproducibility(mmseqs_installed):
    """Test that clustering produces identical results with same random state."""
    import pandas as pd

    from protclust import cluster

    # Create synthetic test data with sequences that are very similar
    # but might get ordered differently without proper controls
    test_sequences = [
        # These sequences differ by just 1-2 amino acids around the similarity threshold
        "MAKVSTAAALLVLAGLVPAAFAADTRTVLVKDGKEIAGFMKGEVAIFNNNGTVEHDAETLCRNAELPDNPISQINPNSIGIF",
        "MAKVSTAAVLLVLAGLVPAAFAADTRTVLVKDGKEIAGFMKGEVAIFNNNGTVEHDAETLCRNAELPDNPISQINPNSIGIF",
        "MAKVSTAAALLVLAGLVPAAFAADTRTVLVKDGKEIAGFMKGEVAIFNNNGTVEHDAETLCRNAELPDNPISQINPNPIGIF",
        # These sequences are shorter but similar - could cause sorting issues
        "MAKVSTAAALLVLAGLVPAAFAADTRTVLVKDGKEIAGFMKGE",
        "MAKVSTAAALLVLAGLVPAAFAADTRTVLVKDGKEIAGFMKGX",  # Non-standard AA to test edge case
        # These have identical length but vary in composition
        "MSTAAALLVLAGLVPAAFAADTRTVLVKDGKEIA",
        "MSTVAVLLVLAGLVPAAFAADTSRVLVKDGKEIA",
        # This one has different capitalization to test case sensitivity
        "mstvavllvlaglvpaafaadtsrvlvkdgkeiaG",
    ]

    test_df = pd.DataFrame(
        {
            "id": [f"seq_{i}" for i in range(len(test_sequences))],
            "sequence": test_sequences,
            "property": [
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
            ],  # Add a property to make sure it's preserved
        }
    )

    # Run clustering twice with the same random state
    result1 = cluster(
        test_df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.9,  # High threshold to make clustering sensitive to small differences
        coverage=0.9,
        random_state=42,  # Fixed random state
    )

    result2 = cluster(
        test_df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.9,
        coverage=0.9,
        random_state=42,  # Same random state
    )

    # Verify exact same cluster assignments for all sequences
    pd.testing.assert_series_equal(
        result1["cluster_representative"], result2["cluster_representative"], check_names=False
    )

    # Count clusters and verify reasonable results
    n_clusters = result1["cluster_representative"].nunique()
    assert 2 <= n_clusters <= 6, f"Expected 2-6 clusters, got {n_clusters}"

    # Optional: Run with different parameters to verify they affect results
    result3 = cluster(
        test_df,
        sequence_col="sequence",
        id_col="id",
        min_seq_id=0.8,  # Lower threshold
        coverage=0.9,
        random_state=42,
    )

    # Verify different parameters produce different clustering
    assert not result1["cluster_representative"].equals(result3["cluster_representative"]), (
        "Changing min_seq_id didn't affect clustering results"
    )

    # Ensure disabling reproducibility controls leads to potential differences
    # We'll run this only if --shuffle 0 is actually implemented in your version
    # This runs the function 5 times WITHOUT controlling randomness to detect potential non-determinism
    try:
        # Monkey patch the check_random_state function to always return None
        import unittest.mock as mock

        import numpy as np

        with mock.patch(
            "protclust.clustering.check_random_state", return_value=np.random.RandomState()
        ):
            results_without_controls = []
            for _ in range(5):
                # Note: If your mmseqs2 function doesn't implement --shuffle 0 yet,
                # this might not catch reproducibility issues
                result = cluster(
                    test_df,
                    sequence_col="sequence",
                    id_col="id",
                    min_seq_id=0.9,
                    coverage=0.9,
                )
                results_without_controls.append(result["cluster_representative"].copy())

            # Check if any of the clustering results differ
            all_same = all(
                results_without_controls[0].equals(r) for r in results_without_controls[1:]
            )

            # If all results are identical despite our attempts to introduce randomness,
            # log a warning but don't fail the test
            if all_same:
                import warnings

                warnings.warn(
                    "All clustering results were identical even without reproducibility controls. "
                    "This test may not effectively detect non-deterministic behavior."
                )
    except Exception as e:
        # If this part fails for any reason, log it but continue
        import warnings

        warnings.warn(f"Optional non-determinism check failed: {str(e)}")
