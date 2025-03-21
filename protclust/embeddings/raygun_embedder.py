"""RayGun embeddings for protein sequences."""

from typing import List, Optional

import numpy as np
import torch

from ..logger import logger
from .baseline import BaseEmbedder


class RayGunEmbedder(BaseEmbedder):
    """Embedder using RayGun model which encodes ESM embeddings."""

    default_pooling = "mean"

    def __init__(
        self,
        esm_model_name: str = "esm2_t33_650M_UR50D",
        device: Optional[str] = None,
    ):
        """
        Initialize the RayGun embedder.

        Args:
            esm_model_name: Name of the ESM model to use for initial embeddings
            device: Device to run inference on ('cpu', 'cuda', 'mps', or None for auto)
        """
        self.esm_model_name = esm_model_name

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
        logger.info(f"Using device: {self.device} for RayGun embeddings")

        # Initialize models
        self._initialize_models()

    def _initialize_models(self):
        """Initialize the ESM and RayGun models."""
        try:
            import esm
        except ImportError:
            raise ImportError("ESM package not found. Install with: pip install fair-esm")

        logger.info(f"Loading ESM model: {self.esm_model_name}")

        # Load ESM model
        self.esm_model, self.esm_alphabet = esm.pretrained.load_model_and_alphabet(
            self.esm_model_name
        )
        self.esm_model = self.esm_model.to(self.device).eval()
        self.batch_converter = self.esm_alphabet.get_batch_converter()

        # Get ESM model layer count for representation extraction
        self.esm_layer = len(self.esm_model.layers)

        logger.info(
            f"ESM model loaded with {sum(p.numel() for p in self.esm_model.parameters()) / 1e6:.1f}M parameters"
        )

        # Load RayGun model
        try:
            logger.info("Loading RayGun model")
            self.ray_model, self.esmtotokdecoder, self.hypparams = torch.hub.load(
                "rohitsinghlab/raygun", "pretrained_uniref50_95000_750M"
            )
            self.ray_model = self.ray_model.to(self.device).eval()

            # Get embedding dimension from the model
            self.embedding_dim = self.ray_model.encoder.output_dim

            logger.info(
                f"RayGun model loaded with {sum(p.numel() for p in self.ray_model.parameters()) / 1e6:.1f}M parameters"
            )
        except Exception as e:
            raise ImportError(f"Failed to load RayGun model: {e}")

    def generate(
        self,
        sequence: str,
        pooling: str = "auto",
        max_length: Optional[int] = 1022,
    ) -> np.ndarray:
        """
        Generate RayGun embedding for a protein sequence.

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

        # Process with ESM model first
        with torch.no_grad():
            # Generate ESM embedding
            batch_labels, batch_strs, batch_tokens = self.batch_converter(
                [("seq", sequence.upper())]
            )
            esm_results = self.esm_model(
                batch_tokens.to(self.device), repr_layers=[self.esm_layer], return_contacts=False
            )
            esm_embedding = esm_results["representations"][self.esm_layer]

            # Remove special tokens (keep only amino acid tokens)
            esm_embedding = esm_embedding[:, 1 : len(sequence) + 1, :]

            # Process with RayGun encoder
            ray_embedding = self.ray_model.encoder(esm_embedding).squeeze().cpu().numpy()

        # RayGun already produces a fixed-length encoding, so we might not need pooling
        # but we'll keep it for consistency with the BaseEmbedder interface
        return self._apply_pooling(ray_embedding, pooling)

    def batch_generate(
        self,
        sequences: List[str],
        pooling: str = "auto",
        max_length: Optional[int] = 1022,
        batch_size: int = 8,
    ) -> List[np.ndarray]:
        """
        Generate embeddings for multiple sequences efficiently.

        Args:
            sequences: List of amino acid sequences
            pooling: Pooling method
            max_length: Maximum sequence length
            batch_size: Number of sequences to process at once

        Returns:
            List of embedding arrays
        """
        if not sequences:
            return []

        results = []

        # Handle empty sequences
        for seq in sequences:
            if not seq:
                if pooling == "none" or (pooling == "auto" and self.default_pooling == "none"):
                    results.append(np.zeros((0, self.embedding_dim)))
                else:
                    results.append(np.zeros(self.embedding_dim))
            else:
                results.append(None)  # Placeholder for non-empty sequences

        # Process in batches
        with torch.no_grad():
            for i in range(0, len(sequences), batch_size):
                batch_seqs = sequences[i : i + batch_size]

                # Process each sequence (truncate if needed)
                processed_batch = []
                for j, seq in enumerate(batch_seqs):
                    if not seq:
                        continue
                    if max_length and len(seq) > max_length:
                        processed_batch.append((f"seq_{j}", seq[:max_length].upper()))
                    else:
                        processed_batch.append((f"seq_{j}", seq.upper()))

                if not processed_batch:
                    continue

                # Generate ESM embeddings for the batch
                batch_labels, batch_strs, batch_tokens = self.batch_converter(processed_batch)
                esm_results = self.esm_model(
                    batch_tokens.to(self.device),
                    repr_layers=[self.esm_layer],
                    return_contacts=False,
                )
                esm_embeddings = esm_results["representations"][self.esm_layer]

                # Process each sequence in the batch
                for j, (label, seq) in enumerate(processed_batch):
                    # Extract ESM embedding for this sequence (remove special tokens)
                    esm_emb = esm_embeddings[j : j + 1, 1 : len(seq) + 1, :]

                    # Process with RayGun encoder
                    ray_emb = self.ray_model.encoder(esm_emb).squeeze().cpu().numpy()

                    # Store result
                    seq_idx = i + int(label.split("_")[1])
                    results[seq_idx] = self._apply_pooling(ray_emb, pooling)

        return results
