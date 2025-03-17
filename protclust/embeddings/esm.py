"""ESM embeddings for protein sequences."""

from typing import List, Optional

import numpy as np
import torch

from ..logger import logger
from .baseline import BaseEmbedder


class ESMEmbedder(BaseEmbedder):
    """Embedder using Facebook/Meta AI's ESM protein language models."""

    default_pooling = "mean"

    def __init__(
        self,
        model_name: str = "esm2_t6_8M_UR50D",
        layer: int = -1,
        device: Optional[str] = None,
    ):
        """
        Initialize the ESM embedder.

        Args:
            model_name: Name of the ESM model to use
            layer: Which layer to extract embeddings from (-1 for last layer)
            device: Device to run inference on ('cpu', 'cuda', 'mps', or None for auto)
        """
        self.model_name = model_name

        # Set device
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
            if (
                device == "cpu"
                and hasattr(torch.backends, "mps")
                and torch.backends.mps.is_available()
            ):
                device = "mps"

        self.device = device
        logger.info(f"Using device: {self.device} for ESM embeddings")

        # Initialize model
        self._initialize_model(layer)

    def _initialize_model(self, layer):
        """Initialize the ESM model and tokenizer."""
        try:
            import esm
        except ImportError:
            raise ImportError("ESM package not found. Install with: pip install fair-esm")

        logger.info(f"Loading ESM model: {self.model_name}")

        self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(self.model_name)
        self.model.to(self.device).eval()

        # Determine the total number of layers in the model
        self.model_layers = len(self.model.layers)

        # Set the layer from which to extract representations
        self.repr_layer = self.model_layers + layer if layer < 0 else layer
        if not 0 <= self.repr_layer < self.model_layers:
            raise ValueError(f"Layer {layer} out of bounds (0-{self.model_layers - 1})")

        # Get embedding dimension for the model
        self.embedding_dim = self.model.embed_tokens.embedding_dim

        # Create the batch converter
        self.batch_converter = self.alphabet.get_batch_converter()

        logger.info(
            f"Model loaded with {sum(p.numel() for p in self.model.parameters()) / 1e6:.1f}M parameters"
        )

    def generate(
        self,
        sequence: str,
        pooling: str = "auto",
        max_length: Optional[int] = 1022,
    ) -> np.ndarray:
        """
        Generate ESM embedding for a protein sequence.

        Args:
            sequence: Amino acid sequence
            pooling: Pooling method
            max_length: Maximum sequence length

        Returns:
            Embedding array
        """
        # Handle empty sequence
        if not sequence:
            if pooling == "none" or (pooling == "auto" and self.default_pooling == "none"):
                return np.zeros((0, self.embedding_dim))
            else:
                return np.zeros(self.embedding_dim)

        # Enforce length limit
        if max_length and len(sequence) > max_length:
            sequence = sequence[:max_length]
            logger.warning(f"Sequence truncated to {max_length} residues")

        # Prepare batch
        batch_data = [("protein", sequence)]

        # Compute embeddings
        with torch.no_grad():
            _, _, batch_tokens = self.batch_converter(batch_data)
            results = self.model(batch_tokens.to(self.device), repr_layers=[self.repr_layer])
            embedding = (
                results["representations"][self.repr_layer][0, 1 : len(sequence) + 1].cpu().numpy()
            )

        return self._apply_pooling(embedding, pooling)

    def batch_generate(
        self,
        sequences: List[str],
        pooling: str = "auto",
        max_length: Optional[int] = 1022,
    ) -> List[np.ndarray]:
        """
        Generate embeddings for multiple sequences efficiently.

        Args:
            sequences: List of amino acid sequences
            pooling: Pooling method
            max_length: Maximum sequence length

        Returns:
            List of embedding arrays
        """
        if not sequences:
            return []

        # Pre-process sequences
        processed_sequences = []
        for i, seq in enumerate(sequences):
            if not seq:
                continue
            if max_length and len(seq) > max_length:
                sequences[i] = seq[:max_length]
            processed_sequences.append((f"protein_{i}", sequences[i]))

        results = []

        # Handle empty sequences
        for seq in sequences:
            if not seq:
                if pooling == "none" or (pooling == "auto" and self.default_pooling == "none"):
                    results.append(np.zeros((0, self.embedding_dim)))
                else:
                    results.append(np.zeros(self.embedding_dim))
            else:
                results.append(None)

        # Batch size
        batch_size = 8

        with torch.no_grad():
            for i in range(0, len(processed_sequences), batch_size):
                batch = processed_sequences[i : i + batch_size]

                _, _, batch_tokens = self.batch_converter(batch)
                outputs = self.model(batch_tokens.to(self.device), repr_layers=[self.repr_layer])
                embeddings = outputs["representations"][self.repr_layer].cpu().numpy()

                for j, (label, seq) in enumerate(batch):
                    seq_idx = int(label.split("_")[1])
                    emb = embeddings[j, 1 : len(seq) + 1]
                    results[seq_idx] = self._apply_pooling(emb, pooling)

        return results
