"""
Data storage for embedding matrices and property scales.
"""

# BLOSUM62 substitution matrix values
# Source: https://www.ncbi.nlm.nih.gov/Class/BLAST/BLOSUM62.txt
BLOSUM62 = {
    "A": [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    "R": [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    "N": [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    "D": [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    "C": [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    "Q": [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    "E": [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    "G": [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    "H": [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    "I": [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    "L": [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    "K": [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    "M": [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    "F": [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    "P": [
        -1,
        -2,
        -2,
        -1,
        -3,
        -1,
        -1,
        -2,
        -2,
        -3,
        -3,
        -1,
        -2,
        -4,
        7,
        -1,
        -1,
        -4,
        -3,
        -2,
    ],
    "S": [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    "T": [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    "W": [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    "Y": [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    "V": [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
    "X": [
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
    ],
}

# BLOSUM90 substitution matrix values
# Source: https://www.ncbi.nlm.nih.gov/Class/BLAST/BLOSUM90.txt
BLOSUM90 = {
    "A": [5, -2, -2, -3, -1, -1, -1, 1, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1],
    "R": [-2, 7, -2, -4, -5, 1, -2, -3, 0, -4, -3, 3, -3, -4, -3, -1, -3, -4, -3, -3],
    "N": [-2, -2, 7, 2, -4, 0, 0, -1, 1, -3, -4, 0, -3, -4, -2, 1, 0, -4, -3, -3],
    "D": [-3, -4, 2, 8, -6, -1, 2, -2, -1, -4, -5, -1, -5, -5, -2, 0, -2, -7, -4, -4],
    "C": [
        -1,
        -5,
        -4,
        -6,
        9,
        -5,
        -7,
        -4,
        -5,
        -2,
        -3,
        -5,
        -3,
        -4,
        -4,
        -2,
        -2,
        -5,
        -4,
        -2,
    ],
    "Q": [-1, 1, 0, -1, -5, 7, 2, -2, 1, -3, -2, 2, -1, -3, -2, 0, -1, -3, -2, -3],
    "E": [-1, -2, 0, 2, -7, 2, 7, -3, -1, -4, -4, -1, -3, -4, -2, 0, -1, -5, -3, -3],
    "G": [1, -3, -1, -2, -4, -2, -3, 7, -3, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -3],
    "H": [-2, 0, 1, -1, -5, 1, -1, -3, 8, -3, -3, -1, -2, -1, -2, -1, -2, -3, 2, -3],
    "I": [-2, -4, -3, -4, -2, -3, -4, -4, -3, 5, 2, -3, 2, 1, -3, -2, -1, -3, -1, 4],
    "L": [-2, -3, -4, -5, -3, -2, -4, -4, -3, 2, 5, -3, 3, 1, -3, -3, -2, -3, -1, 1],
    "K": [-1, 3, 0, -1, -5, 2, -1, -2, -1, -3, -3, 6, -2, -4, -2, 0, -1, -4, -3, -3],
    "M": [-2, -3, -3, -5, -3, -1, -3, -3, -2, 2, 3, -2, 7, 1, -3, -2, -1, -2, -1, 1],
    "F": [-3, -4, -4, -5, -4, -3, -4, -4, -1, 1, 1, -4, 1, 8, -4, -3, -3, 1, 4, -1],
    "P": [
        -1,
        -3,
        -2,
        -2,
        -4,
        -2,
        -2,
        -2,
        -2,
        -3,
        -3,
        -2,
        -3,
        -4,
        10,
        -1,
        -1,
        -4,
        -3,
        -3,
    ],
    "S": [1, -1, 1, 0, -2, 0, 0, 0, -1, -2, -3, 0, -2, -3, -1, 5, 2, -3, -2, -2],
    "T": [0, -3, 0, -2, -2, -1, -1, -2, -2, -1, -2, -1, -1, -3, -1, 2, 5, -3, -2, 0],
    "W": [-4, -4, -4, -7, -5, -3, -5, -3, -3, -3, -3, -4, -2, 1, -4, -3, -3, 15, 2, -4],
    "Y": [-3, -3, -3, -4, -4, -2, -3, -3, 2, -1, -1, -3, -1, 4, -3, -2, -2, 2, 8, -2],
    "V": [-1, -3, -3, -4, -2, -3, -3, -3, -3, 4, 1, -3, 1, -1, -3, -2, 0, -4, -2, 5],
    "X": [
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
    ],
}

# Amino acid property scales
PROPERTY_SCALES = {
    "hydrophobicity": {
        # Kyte & Doolittle (1982)
        # Source: https://doi.org/10.1016/0022-2836(82)90301-1
        "A": 1.8,
        "R": -4.5,
        "N": -3.5,
        "D": -3.5,
        "C": 2.5,
        "Q": -3.5,
        "E": -3.5,
        "G": -0.4,
        "H": -3.2,
        "I": 4.5,
        "L": 3.8,
        "K": -3.9,
        "M": 1.9,
        "F": 2.8,
        "P": -1.6,
        "S": -0.8,
        "T": -0.7,
        "W": -0.9,
        "Y": -1.3,
        "V": 4.2,
        "X": 0.0,
    },
    "charge": {
        # Standard net charge at neutral pH
        "A": 0,
        "R": 1,
        "N": 0,
        "D": -1,
        "C": 0,
        "Q": 0,
        "E": -1,
        "G": 0,
        "H": 0.1,
        "I": 0,
        "L": 0,
        "K": 1,
        "M": 0,
        "F": 0,
        "P": 0,
        "S": 0,
        "T": 0,
        "W": 0,
        "Y": 0,
        "V": 0,
        "X": 0,
    },
    "volume": {
        # Normalized van der Waals volumes (Zamyatnin, 1972)
        # Source: https://doi.org/10.1007/BF01010138
        "A": 0.31,
        "R": 0.85,
        "N": 0.48,
        "D": 0.45,
        "C": 0.45,
        "Q": 0.58,
        "E": 0.56,
        "G": 0.00,
        "H": 0.66,
        "I": 0.76,
        "L": 0.76,
        "K": 0.73,
        "M": 0.76,
        "F": 0.88,
        "P": 0.42,
        "S": 0.31,
        "T": 0.45,
        "W": 1.00,
        "Y": 0.88,
        "V": 0.66,
        "X": 0.50,
    },
    "polarity": {
        # Normalized polarity values (AAindex, Kawashima & Kanehisa, 2000)
        # Source: https://www.genome.jp/aaindex/
        "A": 0.06,
        "R": 1.00,
        "N": 0.65,
        "D": 0.70,
        "C": 0.29,
        "Q": 0.65,
        "E": 0.70,
        "G": 0.11,
        "H": 0.77,
        "I": 0.00,
        "L": 0.00,
        "K": 0.96,
        "M": 0.06,
        "F": 0.10,
        "P": 0.22,
        "S": 0.34,
        "T": 0.34,
        "W": 0.27,
        "Y": 0.41,
        "V": 0.00,
        "X": 0.40,
    },
    "helix_propensity": {
        # Chou & Fasman (1974)
        # Source: https://doi.org/10.1016/S0022-2836(74)80005-7
        "A": 1.42,
        "R": 0.98,
        "N": 0.67,
        "D": 1.01,
        "C": 0.70,
        "Q": 1.11,
        "E": 1.51,
        "G": 0.57,
        "H": 1.00,
        "I": 1.08,
        "L": 1.21,
        "K": 1.14,
        "M": 1.45,
        "F": 1.13,
        "P": 0.57,
        "S": 0.77,
        "T": 0.83,
        "W": 1.08,
        "Y": 0.69,
        "V": 1.06,
        "X": 1.00,
    },
}
