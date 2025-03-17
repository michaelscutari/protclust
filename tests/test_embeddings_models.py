import numpy as np
import pytest


def test_real_esm_embedder():
    """Smoke test with the real ESM model."""
    try:
        from protclust.embeddings import ESMEmbedder
    except ImportError:
        pytest.skip("ESM package not installed")

    embedder = ESMEmbedder(model_name="esm2_t6_8M_UR50D")  # Smallest ESM2 model
    protein_seq = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    embedding = embedder.generate(protein_seq)

    # Verify we get expected shape output
    assert embedding.shape[-1] == embedder.embedding_dim
    assert not np.isnan(embedding).any(), "Embedding contains NaN values"

    # Try batch generation as well
    batch_embeddings = embedder.batch_generate([protein_seq, protein_seq[:20]])
    assert len(batch_embeddings) == 2
    assert batch_embeddings[0].shape[-1] == embedder.embedding_dim
    assert batch_embeddings[1].shape[-1] == embedder.embedding_dim


@pytest.mark.skip("Takes a while - run manually")
def test_real_prottrans_embedder():
    """Smoke test with the real ProtTrans model."""
    try:
        from protclust.embeddings import ProtTransEmbedder
    except ImportError:
        pytest.skip("Transformers package not installed")

    embedder = ProtTransEmbedder(model_name="prot_bert_bfd")  # Or use a smaller model if available
    protein_seq = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    embedding = embedder.generate(protein_seq)

    # Verify we get expected shape output
    assert embedding.shape[-1] == embedder.embedding_dim
    assert not np.isnan(embedding).any(), "Embedding contains NaN values"

    # Try batch generation as well
    batch_embeddings = embedder.batch_generate([protein_seq, protein_seq[:20]])
    assert len(batch_embeddings) == 2
    assert batch_embeddings[0].shape[-1] == embedder.embedding_dim
    assert batch_embeddings[1].shape[-1] == embedder.embedding_dim


@pytest.mark.skip("Requires API access - run manually")
def test_real_esm_api_embedder():
    """Smoke test with the real ESM API."""
    try:
        from protclust.embeddings import ESMAPIEmbedder
    except ImportError:
        pytest.skip("Required packages not installed")

    # Note: This test requires an API key to be set
    try:
        embedder = ESMAPIEmbedder()
        protein_seq = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        embedding = embedder.generate(protein_seq)

        # Verify we get expected output shape
        assert embedding.size > 0
        assert not np.isnan(embedding).any(), "Embedding contains NaN values"
    except Exception as e:
        if "API key" in str(e):
            pytest.skip("ESM API key not configured")
        else:
            raise


def test_real_reduction():
    """Smoke test for dimension reduction with real data."""
    try:
        import numpy as np

        from protclust.embeddings.reduction import reduce_dimensions
    except ImportError:
        pytest.skip("Required packages not installed")

    # Create some random high-dimensional data
    data = np.random.random((100, 1024))

    # Test PCA reduction
    reduced_data, reducer = reduce_dimensions(data, method="pca", n_components=50)

    # Verify output shape and properties
    assert reduced_data.shape == (100, 50)
    assert not np.isnan(reduced_data).any(), "Reduced data contains NaN values"

    # Test applying the reducer to new data
    from protclust.embeddings.reduction import apply_reducer

    new_data = np.random.random((5, 1024))
    transformed = apply_reducer(new_data, reducer)
    assert transformed.shape == (5, 50)
