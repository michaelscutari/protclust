"""Tests for MILP-based splitting functionality with detailed property balance verification."""

from protclust import cluster, milp_split


def test_milp_numeric_properties(realistic_protein_data, mmseqs_installed):
    """Test MILP balancing of multiple numeric properties."""
    df = realistic_protein_data.copy()

    # Cluster the data
    clustered_df = cluster(df, sequence_col="sequence", id_col="id", min_seq_id=0.8)

    # Define numeric properties to balance
    num_props = ["molecular_weight", "hydrophobicity", "isoelectric_point"]

    # Run MILP split balancing multiple numeric properties
    train_df, test_df = milp_split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        balance_cols=num_props,
        balance_weight=1.0,
        time_limit=15,  # Short time limit to keep tests fast
        random_state=42,
    )

    # Verify basic split integrity
    assert len(train_df) + len(test_df) == len(clustered_df)
    assert set(train_df["cluster_representative"]).isdisjoint(
        set(test_df["cluster_representative"])
    )

    # Run a baseline split for comparison
    from protclust import split

    baseline_train, baseline_test = split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        random_state=42,
    )

    # Calculate and compare balance metrics
    milp_imbalances = []
    baseline_imbalances = []

    for prop in num_props:
        # MILP split balance
        milp_train_mean = train_df[prop].mean()
        milp_test_mean = test_df[prop].mean()
        overall_mean = clustered_df[prop].mean()
        milp_imbalance = abs(milp_train_mean - milp_test_mean) / overall_mean
        milp_imbalances.append(milp_imbalance)

        # Baseline split balance
        baseline_train_mean = baseline_train[prop].mean()
        baseline_test_mean = baseline_test[prop].mean()
        baseline_imbalance = abs(baseline_train_mean - baseline_test_mean) / overall_mean
        baseline_imbalances.append(baseline_imbalance)

        # MILP should achieve reasonable balance
        assert milp_imbalance < 0.25, (
            f"Property {prop} poorly balanced: diff ratio = {milp_imbalance:.2f}"
        )

    # MILP should be better than baseline for at least one property
    assert any(m < b for m, b in zip(milp_imbalances, baseline_imbalances)), (
        "MILP split didn't improve balance for any property compared to baseline"
    )


def test_milp_time_limit(realistic_protein_data, mmseqs_installed):
    """Test MILP solver with different time limits."""
    df = realistic_protein_data.copy()

    # Cluster the data
    clustered_df = cluster(df, sequence_col="sequence", id_col="id", min_seq_id=0.8)

    # Run with very short time limit
    import time

    start_time = time.time()
    train_short, test_short = milp_split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        balance_cols=["molecular_weight", "hydrophobicity"],
        time_limit=1,  # Very short - might not find optimal
        random_state=42,
    )
    short_time = time.time() - start_time

    # Run with longer time limit
    start_time = time.time()
    train_long, test_long = milp_split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        balance_cols=["molecular_weight", "hydrophobicity"],
        time_limit=5,  # Still short but more time to optimize
        random_state=42,
    )
    long_time = time.time() - start_time

    # Verify both splits are valid
    assert len(train_short) + len(test_short) == len(clustered_df)
    assert len(train_long) + len(test_long) == len(clustered_df)

    # Longer time should lead to better balancing or the same result
    # Calculate imbalance for short run
    short_imb = (
        abs(train_short["molecular_weight"].mean() - test_short["molecular_weight"].mean())
        / clustered_df["molecular_weight"].mean()
    )

    # Calculate imbalance for long run
    long_imb = (
        abs(train_long["molecular_weight"].mean() - test_long["molecular_weight"].mean())
        / clustered_df["molecular_weight"].mean()
    )

    # Longer run should be at least as good (not significantly worse)
    assert long_imb <= short_imb * 1.1, (
        f"Longer optimization time produced worse balance: {long_imb:.3f} vs {short_imb:.3f}"
    )

    # The timing should reflect the limits
    assert short_time < long_time or abs(short_time - long_time) < 0.5, (
        f"Time limits not reflected in runtime: short={short_time:.2f}s, long={long_time:.2f}s"
    )
