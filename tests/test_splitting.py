import logging

import numpy as np
import pytest

from protclust import (
    cluster,
    milp_split,
    split,
)

logger = logging.getLogger(__name__)


def test_split_data(fluorescence_data):
    """Test splitting data without clustering."""
    # First create a test DataFrame with known groups
    df = fluorescence_data.copy()

    # Assign synthetic groups for testing
    np.random.seed(42)
    n_groups = 10
    df["group"] = np.random.randint(0, n_groups, size=len(df))

    # Run splitting with optimal algorithm
    train_df, test_df = split(df, group_col="group", test_size=0.3, algorithm="optimal")

    # Check results
    assert len(train_df) + len(test_df) == len(df)

    # Check that groups are preserved (no group appears in both train and test)
    train_groups = set(train_df["group"])
    test_groups = set(test_df["group"])
    assert len(train_groups.intersection(test_groups)) == 0

    # Check approximate test size
    test_ratio = len(test_df) / len(df)
    assert 0.25 <= test_ratio <= 0.35  # Allow some flexibility due to group constraints

    # Run splitting with greedy algorithm
    train_df_greedy, test_df_greedy = split(
        df, group_col="group", test_size=0.3, algorithm="greedy"
    )

    # Check results for greedy algorithm
    assert len(train_df_greedy) + len(test_df_greedy) == len(df)

    # Check that groups are preserved (no group appears in both train and test)
    train_groups_greedy = set(train_df_greedy["group"])
    test_groups_greedy = set(test_df_greedy["group"])
    assert len(train_groups_greedy.intersection(test_groups_greedy)) == 0

    # Check approximate test size
    test_ratio_greedy = len(test_df_greedy) / len(df)
    assert 0.25 <= test_ratio_greedy <= 0.35  # Allow some flexibility due to group constraints


def test_milp_split(fluorescence_data, mmseqs_installed):
    """Test MILP-based splitting with balanced properties and distribution similarity."""
    try:
        pass
    except ImportError:
        pytest.skip("PuLP not installed, skipping MILP test")

    df = fluorescence_data.copy()

    # First cluster the data
    clustered_df = cluster(df, sequence_col="sequence", id_col="id", min_seq_id=0.99, coverage=0.8)

    # Run MILP split with distribution similarity, balancing fluorescence
    train_df, test_df = milp_split(
        clustered_df,
        group_col="cluster_representative",
        test_size=0.3,
        balance_cols=["fluorescence"],
        balance_weight=2.0,
        time_limit=10,
    )

    # Check that all samples are accounted for
    assert len(train_df) + len(test_df) == len(clustered_df)

    # Check that groups are preserved (no cluster appears in both splits)
    train_groups = set(train_df["cluster_representative"])
    test_groups = set(test_df["cluster_representative"])
    assert len(train_groups.intersection(test_groups)) == 0

    # Check approximate test size
    test_ratio = len(test_df) / len(clustered_df)
    assert 0.25 <= test_ratio <= 0.35  # Allow some flexibility

    # Calculate distribution statistics for each split
    train_mean = train_df["fluorescence"].mean()
    test_mean = test_df["fluorescence"].mean()
    overall_mean = clustered_df["fluorescence"].mean()

    train_var = train_df["fluorescence"].var()
    test_var = test_df["fluorescence"].var()
    overall_var = clustered_df["fluorescence"].var()

    train_min = train_df["fluorescence"].min()
    test_min = test_df["fluorescence"].min()
    overall_min = clustered_df["fluorescence"].min()

    train_max = train_df["fluorescence"].max()
    test_max = test_df["fluorescence"].max()
    overall_max = clustered_df["fluorescence"].max()

    # Log MILP split distribution results
    logger.info("MILP split distribution statistics for fluorescence:")
    logger.info(f"  Mean: train={train_mean:.4f}, test={test_mean:.4f}, overall={overall_mean:.4f}")
    logger.info(
        f"  Variance: train={train_var:.4f}, test={test_var:.4f}, overall={overall_var:.4f}"
    )
    logger.info(
        f"  Range: train=[{train_min:.4f}, {train_max:.4f}], test=[{test_min:.4f}, {test_max:.4f}], overall=[{overall_min:.4f}, {overall_max:.4f}]"
    )

    # For comparison, perform a naive split and check balance
    naive_train, naive_test = split(clustered_df, group_col="cluster_representative", test_size=0.3)

    # Calculate distribution statistics for naive split
    naive_train_mean = naive_train["fluorescence"].mean()
    naive_test_mean = naive_test["fluorescence"].mean()

    naive_train_var = naive_train["fluorescence"].var()
    naive_test_var = naive_test["fluorescence"].var()

    naive_train_min = naive_train["fluorescence"].min()
    naive_test_min = naive_test["fluorescence"].min()

    naive_train_max = naive_train["fluorescence"].max()
    naive_test_max = naive_test["fluorescence"].max()

    logger.info("Naive split distribution statistics for fluorescence:")
    logger.info(f"  Mean: train={naive_train_mean:.4f}, test={naive_test_mean:.4f}")
    logger.info(f"  Variance: train={naive_train_var:.4f}, test={naive_test_var:.4f}")
    logger.info(
        f"  Range: train=[{naive_train_min:.4f}, {naive_train_max:.4f}], test=[{naive_test_min:.4f}, {naive_test_max:.4f}]"
    )

    # Compute imbalance metrics for MILP and naive splits
    # Mean imbalance
    milp_mean_imbalance = abs(train_mean - test_mean) / overall_mean
    naive_mean_imbalance = abs(naive_train_mean - naive_test_mean) / overall_mean

    # Variance imbalance
    milp_var_imbalance = abs(train_var - test_var) / overall_var
    naive_var_imbalance = abs(naive_train_var - naive_test_var) / overall_var

    # Range coverage imbalance (how well each split captures the full range)
    milp_range_coverage = (train_max - train_min + test_max - test_min) / (
        2 * (overall_max - overall_min)
    )
    naive_range_coverage = (naive_train_max - naive_train_min + naive_test_max - naive_test_min) / (
        2 * (overall_max - overall_min)
    )

    logger.info("Imbalance metrics:")
    logger.info(
        f"  Mean imbalance: MILP={milp_mean_imbalance:.2%}, Naive={naive_mean_imbalance:.2%}"
    )
    logger.info(
        f"  Variance imbalance: MILP={milp_var_imbalance:.2%}, Naive={naive_var_imbalance:.2%}"
    )
    logger.info(
        f"  Range coverage: MILP={milp_range_coverage:.2%}, Naive={naive_range_coverage:.2%}"
    )

    # Assert that MILP split has better distribution balance than naive split
    # We expect at least one of these to be true (though all might not be in every case)
    distribution_improvements = [
        milp_mean_imbalance < naive_mean_imbalance,
        milp_var_imbalance < naive_var_imbalance,
        milp_range_coverage > naive_range_coverage,
    ]

    # Assert that MILP improves at least one distribution metric
    assert any(distribution_improvements), (
        "MILP split did not improve any distribution metrics over naive split"
    )

    # Log which metrics were improved
    if milp_mean_imbalance < naive_mean_imbalance:
        logger.info("MILP improved mean balance")
    if milp_var_imbalance < naive_var_imbalance:
        logger.info("MILP improved variance balance")
    if milp_range_coverage > naive_range_coverage:
        logger.info("MILP improved range coverage")


def test_enhanced_milp_split(fluorescence_data, protein_data=None):
    """Test enhanced MILP-based splitting with categorical and residue-level data."""
    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        pytest.skip("Required packages not installed, skipping enhanced MILP test")

    # If no protein data provided, use a small synthetic dataset
    if protein_data is None:
        # Create synthetic protein data
        protein_data = []
        for i in range(100):
            protein_length = np.random.randint(50, 200)
            protein_data.append(
                {
                    "id": f"protein_{i}",
                    "primary": "".join(
                        np.random.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=protein_length)
                    ),
                    "protein_length": protein_length,
                    "class_label": np.random.randint(0, 5),
                    "fold_label": np.random.randint(0, 10),
                    "superfamily_label": np.random.randint(0, 20),
                    "family_label": np.random.randint(0, 50),
                    # Synthetic residue-level features
                    "secondary_structure": [
                        [np.random.random(), np.random.random(), np.random.random()]
                        for _ in range(protein_length)
                    ],
                    "solvent_accessibility": [
                        [np.random.random(), np.random.random()] for _ in range(protein_length)
                    ],
                }
            )

        protein_df = pd.DataFrame(protein_data)
    else:
        protein_df = pd.DataFrame(protein_data)

    # Add cluster column if it doesn't exist (for testing purposes)
    if "cluster_representative" not in protein_df.columns:
        # Use a simple clustering based on protein length
        protein_df["cluster_representative"] = protein_df["id"].copy()

    # Run enhanced MILP split with protein data
    logger.info("Running enhanced MILP split with protein data")
    train_protein_df, test_protein_df = milp_split(
        protein_df,
        group_col="cluster_representative",
        test_size=0.3,
        balance_cols=["protein_length"],
        balance_weight=2.0,
        time_limit=10,
        random_state=42,
    )

    # Check that all samples are accounted for
    assert len(train_protein_df) + len(test_protein_df) == len(protein_df)

    # Check that groups are preserved (no cluster appears in both splits)
    train_groups = set(train_protein_df["cluster_representative"])
    test_groups = set(test_protein_df["cluster_representative"])
    assert len(train_groups.intersection(test_groups)) == 0

    # Check approximate test size
    test_ratio = len(test_protein_df) / len(protein_df)
    assert 0.25 <= test_ratio <= 0.35  # Allow some flexibility

    # Now test with fluorescence data if provided
    if fluorescence_data is not None:
        df = fluorescence_data.copy()

        # Add some synthetic categorical columns
        df["sequence_type"] = np.random.choice(["wild_type", "mutant", "synthetic"], size=len(df))
        df["expression_level"] = np.random.choice(["low", "medium", "high"], size=len(df))

        # Add some synthetic residue-level data
        df["residue_hydrophobicity"] = [
            np.random.normal(0, 1, size=len(seq)) for seq in df["sequence"]
        ]

        # Add cluster column if it doesn't exist
        if "cluster_representative" not in df.columns:
            df["cluster_representative"] = df["id"].copy()

        logger.info("Running enhanced MILP split with fluorescence data")
        train_df, test_df = milp_split(
            df,
            group_col="cluster_representative",
            test_size=0.3,
            balance_cols=["fluorescence"],
            balance_weight=2.0,
            time_limit=10,
            random_state=42,
        )

        # Check that all samples are accounted for
        assert len(train_df) + len(test_df) == len(df)

        # Calculate distribution statistics for each split
        train_mean = train_df["fluorescence"].mean()
        test_mean = test_df["fluorescence"].mean()
        overall_mean = df["fluorescence"].mean()

        train_var = train_df["fluorescence"].var()
        test_var = test_df["fluorescence"].var()
        overall_var = df["fluorescence"].var()

        # Compute difference in means as a percentage of overall mean
        mean_diff_pct = abs(train_mean - test_mean) / overall_mean * 100
        var_diff_pct = abs(train_var - test_var) / overall_var * 100

        # Log results
        logger.info("Fluorescence distribution balance:")
        logger.info(
            f"  Mean: train={train_mean:.4f}, test={test_mean:.4f}, diff={mean_diff_pct:.2f}%"
        )
        logger.info(
            f"  Variance: train={train_var:.4f}, test={test_var:.4f}, diff={var_diff_pct:.2f}%"
        )

        # Check categorical balance
        for cat_col in ["sequence_type", "expression_level"]:
            # Compare category distributions
            train_dist = train_df[cat_col].value_counts(normalize=True)
            test_dist = test_df[cat_col].value_counts(normalize=True)

            # Compute Jensen-Shannon divergence (simplified)
            all_cats = set(train_dist.index).union(set(test_dist.index))
            js_div = 0
            for cat in all_cats:
                p = train_dist.get(cat, 0)
                q = test_dist.get(cat, 0)
                js_div += abs(p - q)
            js_div /= 2  # Normalize

            logger.info(f"Categorical balance for {cat_col}:")
            logger.info(f"  Train distribution: {dict(train_dist)}")
            logger.info(f"  Test distribution: {dict(test_dist)}")
            logger.info(f"  JS divergence: {js_div:.4f} (lower is better)")

            # Assert reasonable balance of categories
            assert js_div < 0.5, (
                f"Poor categorical balance for {cat_col}: JS divergence = {js_div:.4f}"
            )

        # For comparison, perform a naive split and check balance
        from random import Random

        rng = Random(42)
        test_indices = set(rng.sample(range(len(df)), int(0.3 * len(df))))
        naive_test = df.iloc[list(test_indices)]
        naive_train = df.iloc[list(set(range(len(df))) - test_indices)]

        naive_train_mean = naive_train["fluorescence"].mean()
        naive_test_mean = naive_test["fluorescence"].mean()
        naive_mean_diff_pct = abs(naive_train_mean - naive_test_mean) / overall_mean * 100

        logger.info(f"Naive split mean difference: {naive_mean_diff_pct:.2f}%")

        # The MILP split should generally achieve better balance than naive random split
        # (though this might not always be true due to clustering constraints)
        logger.info(
            f"MILP vs Naive mean balance improvement: {naive_mean_diff_pct - mean_diff_pct:.2f}%"
        )


def test_empty_dataset_handling():
    """Test splitting functions with empty datasets."""
    import pandas as pd

    from protclust import split

    # Create empty DataFrame with the required columns
    empty_df = pd.DataFrame(columns=["cluster_representative", "id"])

    # Test basic split - should raise ValueError for empty DataFrame
    with pytest.raises(ValueError, match="Cannot split empty DataFrame"):
        split(empty_df, group_col="cluster_representative")


def test_milp_split_categorical_handling(fluorescence_data):
    """Test MILP splitting with categorical variables."""
    import numpy as np

    from protclust import milp_split

    df = fluorescence_data.head(20).copy()

    # Add categorical columns
    df["category"] = np.random.choice(["A", "B", "C"], size=len(df))
    df["binary"] = np.random.choice([0, 1], size=len(df))

    # Add cluster column simulation
    df["cluster_representative"] = df.index.astype(str)

    try:
        # Test with categorical columns
        train_df, test_df = milp_split(
            df,
            group_col="cluster_representative",
            test_size=0.3,
            time_limit=1,  # Short timeout to avoid hanging tests
        )

        # Check that split was created
        assert len(train_df) + len(test_df) == len(df)
    except ImportError:
        # Skip if PuLP not available
        pytest.skip("PuLP not installed, skipping MILP test")
