from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .logger import logger
from .utils import check_random_state


def split(
    df: pd.DataFrame,
    group_col: str = "cluster_representative",
    test_size: float = 0.2,
    random_state: Optional[Union[int, np.random.RandomState]] = None,
    tolerance: float = 0.05,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits DataFrame into train/test sets based on grouping in a specified column.

    This function ensures that all sequences from the same cluster end up in the same split,
    preventing data leakage between train and test sets.

    Parameters:
        df (pd.DataFrame): DataFrame to split (should be pre-clustered).
        group_col (str): Column by which to group before splitting (default: "cluster_representative").
        test_size (float): Desired fraction of data in test set (default 0.2).
        random_state (int or RandomState, optional): Random state for reproducibility.
        tolerance (float): Acceptable deviation from test_size (default 0.05).

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: (train_df, test_df)

    Raises:
        ValueError: If the achieved test size deviates too much from the target.
    """
    logger.info(f"Splitting data by '{group_col}' with target test size {test_size}")

    # Get random state for reproducibility
    rng = check_random_state(random_state)

    if group_col not in df.columns:
        raise ValueError(
            f"Column '{group_col}' not found in DataFrame. Did you run clustering first?"
        )

    total_sequences = len(df)
    if total_sequences == 0:
        raise ValueError("Cannot split empty DataFrame")

    target_test_count = int(round(test_size * total_sequences))

    logger.info(f"Total sequence count: {total_sequences}")
    logger.info(f"Target test count: {target_test_count}")

    # Get group sizes
    size_per_group = df.groupby(group_col).size()

    # Create tuples of (group, size) and sort them deterministically
    group_size_pairs = [(group, size) for group, size in size_per_group.items()]
    group_size_pairs.sort(key=lambda x: str(x[0]))  # Deterministic sorting

    # Shuffle using the random state
    rng.shuffle(group_size_pairs)

    # Use greedy algorithm to select groups for test set
    chosen_groups = _greedy_subset_selection(group_size_pairs, target_test_count)

    logger.debug(f"Selected {len(chosen_groups)} groups for test set")

    test_df = df[df[group_col].isin(chosen_groups)]
    train_df = df[~df[group_col].isin(chosen_groups)]

    achieved_test_fraction = len(test_df) / total_sequences

    logger.info(f"Train set: {len(train_df)} sequences ({len(train_df) / total_sequences:.2%})")
    logger.info(f"Test set: {len(test_df)} sequences ({achieved_test_fraction:.2%})")

    if abs(achieved_test_fraction - test_size) > tolerance:
        raise ValueError(
            f"Could not achieve target test size {test_size:.1%}. "
            f"Achieved {achieved_test_fraction:.1%} (tolerance: ±{tolerance:.1%}). "
            f"Consider adjusting tolerance or clustering parameters."
        )

    return train_df, test_df


def _greedy_subset_selection(group_size_pairs: List[Tuple], target_count: int) -> List:
    """
    Use a greedy algorithm to find subset of groups that gets close to target count.

    Parameters:
        group_size_pairs (list): List of (group, size) tuples.
        target_count (int): Target number of sequences.

    Returns:
        list: List of group identifiers chosen.
    """
    logger.debug(f"Using greedy algorithm for subset selection with {len(group_size_pairs)} groups")

    # Sort by size (descending) for better packing
    sorted_pairs = sorted(group_size_pairs, key=lambda x: x[1], reverse=True)

    selected_groups = []
    current_count = 0

    for group, size in sorted_pairs:
        if current_count + size <= target_count:
            selected_groups.append(group)
            current_count += size
        elif current_count == 0:
            # If we can't fit any group, take the smallest one
            selected_groups.append(group)
            current_count += size
            break

    # If we're still far from target, try to add smaller groups
    if current_count < target_count * 0.8:  # If we're using less than 80% of target
        remaining_pairs = [(g, s) for g, s in group_size_pairs if g not in selected_groups]
        remaining_pairs.sort(key=lambda x: x[1])  # Sort by size (ascending)

        for group, size in remaining_pairs:
            if current_count + size <= target_count * 1.2:  # Don't exceed 120% of target
                selected_groups.append(group)
                current_count += size

    return selected_groups


def milp_split(
    df: pd.DataFrame,
    group_col: str = "cluster_representative",
    test_size: float = 0.2,
    balance_cols: Optional[List[str]] = None,
    random_state: Optional[Union[int, np.random.RandomState]] = None,
    time_limit: int = 60,
    balance_weight: float = 1.0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split DataFrame using Mixed Integer Linear Programming for optimal balancing.

    This function uses MILP to find the optimal split that balances numerical columns
    between train and test sets while respecting cluster boundaries.

    Parameters:
        df (pd.DataFrame): DataFrame to split (should be pre-clustered).
        group_col (str): Column by which to group before splitting.
        test_size (float): Desired fraction of data in test set.
        balance_cols (List[str], optional): Numerical columns to balance between splits.
        random_state (int or RandomState, optional): Random state for reproducibility.
        time_limit (int): Maximum time in seconds for MILP solver.
        balance_weight (float): Weight for balancing terms in objective function.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: (train_df, test_df)

    Note:
        Requires 'pulp' package. Falls back to regular split() if not available.
    """
    try:
        from pulp import (
            PULP_CBC_CMD,
            LpBinary,
            LpMinimize,
            LpProblem,
            LpStatus,
            LpVariable,
            lpSum,
            value,
        )
    except ImportError:
        logger.warning("PuLP not available. Falling back to regular split()")
        return split(df, group_col, test_size, random_state)

    logger.info(f"Performing MILP-based splitting with target test size {test_size}")

    if balance_cols is None:
        balance_cols = []

    # Validate balance columns
    balance_cols = [col for col in balance_cols if col in df.columns]
    if not balance_cols:
        logger.warning("No valid balance columns found. Using regular split.")
        return split(df, group_col, test_size, random_state)

    # Get group information
    group_sizes = df.groupby(group_col).size()
    groups = list(group_sizes.index)
    total_size = len(df)
    target_test_count = test_size * total_size

    logger.info(f"Balancing {len(balance_cols)} columns across {len(groups)} groups")

    # Calculate property sums per group
    group_properties = {}
    property_totals = {}

    for col in balance_cols:
        # Only use numeric columns
        if not pd.api.types.is_numeric_dtype(df[col]):
            logger.warning(f"Column '{col}' is not numeric, skipping")
            continue

        group_sums = df.groupby(group_col)[col].sum()
        group_properties[col] = group_sums
        property_totals[col] = df[col].sum()

    if not group_properties:
        logger.warning("No numeric balance columns found. Using regular split.")
        return split(df, group_col, test_size, random_state)

    # Create MILP problem
    prob = LpProblem("ClusterSplit", LpMinimize)

    # Decision variables: 1 if group is in test set, 0 if in train set
    group_vars = {group: LpVariable(f"group_{i}", cat=LpBinary) for i, group in enumerate(groups)}

    # Variables for absolute deviations
    size_deviation = LpVariable("size_deviation", lowBound=0)
    balance_deviations = {
        col: LpVariable(f"balance_{col}", lowBound=0) for col in group_properties.keys()
    }

    # Objective: minimize size deviation + weighted balance deviations
    objective_terms = [size_deviation]
    if balance_deviations:
        objective_terms.append(balance_weight * lpSum(balance_deviations.values()))

    prob += lpSum(objective_terms)

    # Constraint: size deviation
    test_count = lpSum([group_sizes[group] * group_vars[group] for group in groups])
    prob += size_deviation >= test_count - target_test_count
    prob += size_deviation >= target_test_count - test_count

    # Constraints: property balance
    for col, group_prop in group_properties.items():
        test_prop_sum = lpSum([group_vars[group] * group_prop[group] for group in groups])
        target_test_prop = test_size * property_totals[col]

        prob += balance_deviations[col] >= test_prop_sum - target_test_prop
        prob += balance_deviations[col] >= target_test_prop - test_prop_sum

    # Solve the problem
    try:
        solver = PULP_CBC_CMD(timeLimit=time_limit, msg=False)
        prob.solve(solver)
    except Exception as e:
        logger.warning(f"MILP solver failed: {e}. Falling back to regular split.")
        return split(df, group_col, test_size, random_state)

    logger.info(f"MILP solution status: {LpStatus[prob.status]}")

    if prob.status != 1:  # Not optimal
        logger.warning("MILP solution not optimal. Falling back to regular split.")
        return split(df, group_col, test_size, random_state)

    # Extract solution
    test_groups = [group for group, var in group_vars.items() if value(var) > 0.5]

    # Create splits
    test_df = df[df[group_col].isin(test_groups)]
    train_df = df[~df[group_col].isin(test_groups)]

    # Report results
    achieved_test_fraction = len(test_df) / total_size
    logger.info(f"Train set: {len(train_df)} sequences ({len(train_df) / total_size:.2%})")
    logger.info(f"Test set: {len(test_df)} sequences ({achieved_test_fraction:.2%})")

    # Report balance for each column
    for col in group_properties.keys():
        train_mean = train_df[col].mean() if len(train_df) > 0 else 0
        test_mean = test_df[col].mean() if len(test_df) > 0 else 0
        diff_pct = abs(test_mean - train_mean) / (abs(train_mean) if train_mean != 0 else 1) * 100
        logger.info(
            f"Column '{col}': train_mean={train_mean:.3f}, test_mean={test_mean:.3f}, diff={diff_pct:.1f}%"
        )

    return train_df, test_df
