"""Protein sequence embedding functionality for MMseqsPy."""

from .api import (
    embed_sequences,
    get_embeddings,
    list_available_embedders,
    register_embedder,
)
from .baseline import (
    AACompositionEmbedder,
    BLOSUM90Embedder,
    BLOSUMEmbedder,
    DiAACompositionEmbedder,
    OneHotEmbedder,
    PropertyEmbedder,
)
from .esm import ESMEmbedder
from .reduction import apply_reducer, load_reducer, reduce_dimensions, save_reducer
from .remote import ESMAPIEmbedder
from .storage import (
    get_embeddings_from_df,
    get_embeddings_from_hdf,
    list_embeddings_in_hdf,
    store_embeddings_in_df,
    store_embeddings_in_hdf,
)
from .transformers import ProtTransEmbedder

# Register embedders
register_embedder("esm", ESMEmbedder)
register_embedder("prottrans", ProtTransEmbedder)
register_embedder("esm_api", ESMAPIEmbedder)


# Convenience functions for common embedding types
def blosum62(df, sequence_col="sequence", **kwargs):
    """Add BLOSUM62 embeddings to DataFrame."""
    return embed_sequences(df, "blosum62", sequence_col, **kwargs)


def blosum90(df, sequence_col="sequence", **kwargs):
    """Add BLOSUM90 embeddings to DataFrame."""
    return embed_sequences(df, "blosum90", sequence_col, **kwargs)


def aac(df, sequence_col="sequence", **kwargs):
    """Add amino acid composition embeddings to DataFrame."""
    return embed_sequences(df, "aac", sequence_col, **kwargs)


def property_embedding(df, sequence_col="sequence", **kwargs):
    """Add amino acid property embeddings to DataFrame."""
    return embed_sequences(df, "property", sequence_col, **kwargs)


def onehot(df, sequence_col="sequence", **kwargs):
    """Add one-hot encoded embeddings to DataFrame."""
    return embed_sequences(df, "onehot", sequence_col, **kwargs)


def esm2(df, sequence_col="sequence", model_name="esm2_t6_8M_UR50D", **kwargs):
    """Add ESM embeddings to DataFrame."""
    return embed_sequences(
        df, "esm", sequence_col, model_kwargs={"model_name": model_name, **kwargs}
    )


def prot_bert(df, sequence_col="sequence", **kwargs):
    """Add ProtBERT embeddings to DataFrame."""
    return embed_sequences(
        df, "prottrans", sequence_col, model_kwargs={"model_name": "bert", **kwargs}
    )


def prot_t5(df, sequence_col="sequence", **kwargs):
    """Add ProtT5 embeddings to DataFrame."""
    return embed_sequences(
        df, "prottrans", sequence_col, model_kwargs={"model_name": "t5", **kwargs}
    )


def esm_api(df, sequence_col="sequence", **kwargs):
    """Add ESM API embeddings to DataFrame."""
    return embed_sequences(df, "esm_api", sequence_col, model_kwargs={**kwargs})


__all__ = [
    "embed_sequences",
    "get_embeddings",
    "list_available_embedders",
    "register_embedder",
    "BLOSUMEmbedder",
    "AACompositionEmbedder",
    "PropertyEmbedder",
    "OneHotEmbedder",
    "BLOSUM90Embedder",
    "DiAACompositionEmbedder",
    "reduce_dimensions",
    "apply_reducer",
    "save_reducer",
    "load_reducer",
    "store_embeddings_in_df",
    "store_embeddings_in_hdf",
    "get_embeddings_from_df",
    "get_embeddings_from_hdf",
    "list_embeddings_in_hdf",
    "blosum62",
    "blosum90",
    "aac",
    "property_embedding",
    "onehot",
]
