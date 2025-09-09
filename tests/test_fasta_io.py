"""Tests for FASTA file I/O functionality."""

import os
import tempfile

import pandas as pd
import pytest

from protclust.utils import load_fasta_to_dataframe, save_dataframe_to_fasta


def test_load_fasta_to_dataframe():
    """Test loading a FASTA file into a DataFrame."""
    # Create a temporary FASTA file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_file:
        temp_file.write(">seq1\nACGT\n")
        temp_file.write(">seq2\nTGCA\n")
        temp_file.write(">seq3 description\nAGCT\nAGCT\n")  # Multi-line sequence
        temp_path = temp_file.name

    try:
        # Test loading the file
        df = load_fasta_to_dataframe(temp_path)

        # Check the DataFrame structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert list(df.columns) == ["id", "sequence"]

        # Check the values
        assert df.loc[0, "id"] == "seq1"
        assert df.loc[0, "sequence"] == "ACGT"
        assert df.loc[1, "id"] == "seq2"
        assert df.loc[1, "sequence"] == "TGCA"
        assert df.loc[2, "id"] == "seq3 description"
        assert df.loc[2, "sequence"] == "AGCTAGCT"
    finally:
        # Clean up
        os.unlink(temp_path)


def test_load_fasta_file_not_found():
    """Test error handling when file doesn't exist."""
    non_existent_file = "non_existent_file.fasta"
    with pytest.raises(FileNotFoundError):
        load_fasta_to_dataframe(non_existent_file)


def test_load_empty_fasta():
    """Test error handling with an empty FASTA file."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_file:
        temp_path = temp_file.name

    try:
        with pytest.raises(ValueError):
            load_fasta_to_dataframe(temp_path)
    finally:
        os.unlink(temp_path)


def test_save_dataframe_to_fasta():
    """Test saving a DataFrame to a FASTA file."""
    # Create a DataFrame with sequences
    data = {"id": ["seq1", "seq2", "seq3"], "sequence": ["ACGT", "TGCA", "AGCT"]}
    df = pd.DataFrame(data)

    # Create a temporary file path
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        fasta_path = tmp.name

    try:
        # Save the DataFrame to FASTA
        save_dataframe_to_fasta(df, fasta_path)

        # Check if the file exists
        assert os.path.exists(fasta_path)

        # Verify the content
        with open(fasta_path, "r") as f:
            content = f.read()

        expected_content = ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAGCT\n"
        assert content == expected_content

        # Test loading the saved file back
        loaded_df = load_fasta_to_dataframe(fasta_path)
        assert len(loaded_df) == 3
        assert list(loaded_df["id"]) == ["seq1", "seq2", "seq3"]
        assert list(loaded_df["sequence"]) == ["ACGT", "TGCA", "AGCT"]
    finally:
        # Clean up
        os.unlink(fasta_path)


def test_save_with_custom_columns():
    """Test saving with custom column names."""
    # Create a DataFrame with different column names
    data = {"seq_id": ["seq1", "seq2"], "protein_sequence": ["ACGT", "TGCA"]}
    df = pd.DataFrame(data)

    # Create a temporary file path
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        fasta_path = tmp.name

    try:
        # Save with custom column names
        save_dataframe_to_fasta(df, fasta_path, id_col="seq_id", sequence_col="protein_sequence")

        # Verify by loading
        loaded_df = load_fasta_to_dataframe(fasta_path)
        assert len(loaded_df) == 2
        assert list(loaded_df["id"]) == ["seq1", "seq2"]
        assert list(loaded_df["sequence"]) == ["ACGT", "TGCA"]
    finally:
        # Clean up
        os.unlink(fasta_path)


def test_save_empty_dataframe():
    """Test error handling when saving an empty DataFrame."""
    empty_df = pd.DataFrame({"id": [], "sequence": []})

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        fasta_path = tmp.name

    try:
        with pytest.raises(ValueError):
            save_dataframe_to_fasta(empty_df, fasta_path)
    finally:
        os.unlink(fasta_path)


def test_save_with_missing_columns():
    """Test error handling when columns are missing."""
    # Missing sequence column
    df1 = pd.DataFrame({"id": ["seq1"]})

    # Missing ID column
    df2 = pd.DataFrame({"sequence": ["ACGT"]})

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        fasta_path = tmp.name

    try:
        with pytest.raises(ValueError):
            save_dataframe_to_fasta(df1, fasta_path)

        with pytest.raises(ValueError):
            save_dataframe_to_fasta(df2, fasta_path)
    finally:
        os.unlink(fasta_path)
