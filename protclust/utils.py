import os
import shutil

import numpy as np
import pandas as pd

from .logger import logger


def check_random_state(random_state):
    """
    Turn seed into a np.random.RandomState instance.

    Parameters:
        random_state (None, int, or RandomState): If None, return the global
            RandomState. If int, return a new RandomState with the seed.
            If RandomState, return it unchanged.

    Returns:
        RandomState: NumPy RandomState object
    """
    if random_state is None:
        return np.random.RandomState()
    elif isinstance(random_state, (int, np.integer)):
        return np.random.RandomState(random_state)
    elif isinstance(random_state, np.random.RandomState):
        return random_state
    raise ValueError(f"random_state must be None, int, or RandomState, got {type(random_state)}")


def load_fasta_to_dataframe(fasta_path):
    """
    Load sequences from a FASTA file into a pandas DataFrame.

    Parameters:
        fasta_path (str): Path to the FASTA file.

    Returns:
        pd.DataFrame: DataFrame with 'id' and 'sequence' columns.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the FASTA file is empty or invalid.
    """
    if not os.path.exists(fasta_path):
        logger.error(f"FASTA file not found: {fasta_path}")
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    ids = []
    sequences = []
    current_id = None
    current_seq = []

    logger.debug(f"Loading sequences from FASTA file: {fasta_path}")

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous sequence if there was one
                if current_id is not None:
                    ids.append(current_id)
                    sequences.append("".join(current_seq))

                # Start a new sequence
                current_id = line[1:].strip()  # Remove the '>' character
                current_seq = []
            else:
                # Add to the current sequence
                current_seq.append(line)

    # Don't forget to save the last sequence
    if current_id is not None:
        ids.append(current_id)
        sequences.append("".join(current_seq))

    if not ids:
        logger.error(f"No sequences found in FASTA file: {fasta_path}")
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

    logger.info(f"Loaded {len(ids)} sequences from FASTA file")

    return pd.DataFrame({"id": ids, "sequence": sequences})


def save_dataframe_to_fasta(df, fasta_path, id_col="id", sequence_col="sequence"):
    """
    Save sequences from a DataFrame to a FASTA file.

    Parameters:
        df (pd.DataFrame): DataFrame containing the sequences.
        fasta_path (str): Path where the FASTA file will be saved.
        id_col (str): Name of the column containing sequence IDs. Default is 'id'.
        sequence_col (str): Name of the column containing sequences. Default is 'sequence'.

    Raises:
        ValueError: If required columns are missing or if the DataFrame is empty.
    """
    if df.empty:
        logger.error("Cannot save empty DataFrame to FASTA file")
        raise ValueError("DataFrame is empty")

    if id_col not in df.columns:
        logger.error(f"ID column '{id_col}' not found in DataFrame")
        raise ValueError(f"ID column '{id_col}' not found in DataFrame")

    if sequence_col not in df.columns:
        logger.error(f"Sequence column '{sequence_col}' not found in DataFrame")
        raise ValueError(f"Sequence column '{sequence_col}' not found in DataFrame")

    logger.debug(f"Saving {len(df)} sequences to FASTA file: {fasta_path}")

    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(fasta_path), exist_ok=True)

    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            f.write(f">{row[id_col]}\n{row[sequence_col]}\n")

    logger.info(f"Saved {len(df)} sequences to FASTA file: {fasta_path}")


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
