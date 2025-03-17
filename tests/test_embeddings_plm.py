"""Tests for protein embedders."""

import json
import sys
import tempfile
from pathlib import Path

import pytest
import torch


# Set up mocks for external dependencies
@pytest.fixture(autouse=True)
def mock_esm(monkeypatch):
    """Mock ESM library."""

    class MockModel:
        def __init__(self):
            self.num_layers = 5
            self.layers = [None] * self.num_layers  # Add this line
            self.args = type("obj", (object,), {"embed_dim": 320})
            self.embed_tokens = type("obj", (object,), {"embedding_dim": 320})()

        def __call__(self, tokens, repr_layers):
            batch_size, seq_len = tokens.shape
            representations = {
                layer: torch.ones((batch_size, seq_len, 320)) for layer in repr_layers
            }
            return {"representations": representations}

        def eval(self):
            return self

        def to(self, device):
            return self

        def parameters(self):
            class MockParameter:
                def numel(self):
                    return 1000000

            return [MockParameter()]

    class MockAlphabet:
        def get_batch_converter(self):
            def converter(batch_data):
                batch_labels = []
                batch_strs = []
                # Calculate max sequence length and add 2 for BOS and EOS tokens
                max_len = max(len(seq) for _, seq in batch_data) + 2
                batch_tokens = torch.zeros((len(batch_data), max_len), dtype=torch.long)

                for i, (label, seq) in enumerate(batch_data):
                    batch_labels.append(label)
                    batch_strs.append(seq)
                    batch_tokens[i, 1 : len(seq) + 1] = torch.ones(len(seq))

                return batch_labels, batch_strs, batch_tokens

            return converter

    # Patch esm module
    mock_module = type(
        "MockESM",
        (),
        {
            "pretrained": type(
                "MockPretrained",
                (),
                {
                    "load_model_and_alphabet": lambda model_name: (
                        MockModel(),
                        MockAlphabet(),
                    )
                },
            )
        },
    )

    monkeypatch.setattr("sys.modules", {**sys.modules, "esm": mock_module})

    return mock_module


@pytest.fixture(autouse=True)
def mock_transformers(monkeypatch):
    """Mock transformers library."""

    class MockConfig:
        def __init__(self, model_type="bert"):
            self.hidden_size = 16
            self.num_hidden_layers = 5
            self.num_layers = 5

    class MockT5Encoder:
        def __init__(self):
            self.config = MockConfig("t5")

    class MockBertModel:
        def __init__(self, **kwargs):
            self.config = MockConfig()

        def __call__(self, **kwargs):
            batch_size = kwargs["input_ids"].shape[0]
            seq_len = kwargs["input_ids"].shape[1]
            hidden_states = tuple(torch.ones((batch_size, seq_len, 16)) for _ in range(6))
            return type(
                "obj",
                (object,),
                {
                    "hidden_states": hidden_states,
                    "last_hidden_state": hidden_states[-1],
                },
            )

        def eval(self):
            return self

        def to(self, device):
            return self

        def parameters(self):
            # Return an iterable with mock parameters
            class MockParameter:
                def numel(self):
                    return 1000000  # 1M params for simplicity

            return [MockParameter()]

    class MockT5Model:
        def __init__(self, **kwargs):
            self.config = MockConfig("t5")
            self.encoder = MockT5Encoder()

        def __call__(self, **kwargs):
            batch_size = kwargs["input_ids"].shape[0]
            seq_len = kwargs["input_ids"].shape[1]
            last_hidden_state = torch.ones((batch_size, seq_len, 16))
            return type("obj", (object,), {"last_hidden_state": last_hidden_state})

        def eval(self):
            return self

        def to(self, device):
            return self

        def parameters(self):
            # Return an iterable with mock parameters
            class MockParameter:
                def numel(self):
                    return 1000000  # 1M params for simplicity

            return [MockParameter()]

    class MockTokenizer:
        def __init__(self, *args, **kwargs):
            pass

        def __call__(self, texts, **kwargs):
            if isinstance(texts, str):
                texts = [texts]
            batch_size = len(texts)
            seq_len = max(len(text.split()) for text in texts) + 2
            return {
                "input_ids": torch.ones((batch_size, seq_len), dtype=torch.long),
                "attention_mask": torch.ones((batch_size, seq_len), dtype=torch.long),
            }

    # Patch transformers module
    mock_module = type(
        "MockTransformers",
        (),
        {
            "AutoTokenizer": type(
                "MockAutoTokenizer",
                (),
                {"from_pretrained": lambda *args, **kwargs: MockTokenizer()},
            ),
            "AutoModel": type(
                "MockAutoModel",
                (),
                {"from_pretrained": lambda *args, **kwargs: MockBertModel()},
            ),
            "T5EncoderModel": type(
                "MockT5EncoderModel",
                (),
                {"from_pretrained": lambda *args, **kwargs: MockT5Model()},
            ),
        },
    )

    monkeypatch.setattr("sys.modules", {**sys.modules, "transformers": mock_module})

    return mock_module


@pytest.fixture
def mock_esm_api(monkeypatch):
    class MockResponse:
        def __init__(self, status_code=200):
            self.status_code = status_code
            self.text = json.dumps(
                {
                    "representations": {
                        "per_tok": [[0.1] * 32] * 10  # 10 tokens, 32 dimensions
                    }
                }
            )

        def json(self):
            return json.loads(self.text)

        def raise_for_status(self):
            if self.status_code >= 400:
                raise Exception(f"HTTP Error: {self.status_code}")

    # Create a mock post function
    def mock_post(*args, **kwargs):
        return MockResponse()

    # Use autouse=True to apply this fixture automatically to all tests that need it
    # Mock at ESMAPIEmbedder.make_api_request level instead of the requests module
    monkeypatch.setattr(
        "protclust.embeddings.remote.ESMAPIEmbedder._make_request",
        lambda self, url, data: MockResponse(),
    )


@pytest.fixture
def temp_cache_dir():
    """Create temporary directory for caching."""
    with tempfile.TemporaryDirectory() as tempdir:
        yield tempdir


def test_esm_embedder(mock_esm):
    """Test the ESM embedder basic functionality."""
    from protclust.embeddings import ESMEmbedder

    # Initialize embedder
    embedder = ESMEmbedder(model_name="esm2_t6_8M_UR50D")

    # Generate embedding for a sequence
    embedding = embedder.generate("ACDEFGHIKL")  # 320-dim since esm2-8m is default
    assert embedding.shape == (320,)  # Default pooled embedding

    # Test per-residue embeddings
    embedding = embedder.generate("ACDEFGHIKL", pooling="none")
    assert embedding.shape == (10, 320)  # 10 residues, 16 dimensions

    # Test empty sequence
    embedding = embedder.generate("")
    assert embedding.shape == (320,)  # Default pooled for empty

    # Test batch generation
    embeddings = embedder.batch_generate(["ACDEFG", "KLM", ""])
    assert len(embeddings) == 3
    assert embeddings[0].shape == (320,)
    assert embeddings[1].shape == (320,)
    assert embeddings[2].shape == (320,)


def test_prottrans_embedder(mock_transformers):
    """Test the ProtTrans embedder basic functionality."""
    from protclust.embeddings import ProtTransEmbedder

    # Test BERT model
    bert_embedder = ProtTransEmbedder(model_name="bert")

    # Generate embedding
    embedding = bert_embedder.generate("ACDEFGHIKL")
    assert embedding.shape == (16,)  # Default pooled embedding

    # Test T5 model
    t5_embedder = ProtTransEmbedder(model_name="t5")

    # Generate embedding
    embedding = t5_embedder.generate("ACDEFGHIKL")
    assert embedding.shape == (16,)

    # Test batch generation
    embeddings = bert_embedder.batch_generate(["ACDEFG", "KLM", ""])
    assert len(embeddings) == 3
    assert all(emb.shape == (16,) for emb in embeddings)


@pytest.mark.skip("Skipping due to requests/charset_normalizer circular import issue")
def test_esm_api_embedder(mock_esm_api, temp_cache_dir):
    """Test the ESM API embedder basic functionality."""
    from protclust.embeddings.remote import ESMAPIEmbedder

    # Initialize embedder with cache
    embedder = ESMAPIEmbedder(cache_dir=temp_cache_dir)

    # Generate embedding
    embedding = embedder.generate("ACDEFGHIKL")
    assert embedding.shape == (32,)  # From mock response

    # Test caching
    cache_files = list(Path(temp_cache_dir).glob("*.npz"))
    assert len(cache_files) == 1

    # Second call should use cache
    embedding2 = embedder.generate("ACDEFGHIKL")
    assert embedding2.shape == (32,)

    # Test batch generation
    embeddings = embedder.batch_generate(["ACDEFG", "KLM"])
    assert len(embeddings) == 2
    assert all(emb.shape == (32,) for emb in embeddings)
