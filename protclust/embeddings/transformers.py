"""Protein embeddings using transformer models from HuggingFace."""

from typing import List, Optional

import numpy as np
import torch

from ..logger import logger
from .baseline import BaseEmbedder


class ProtTransEmbedder(BaseEmbedder):
    """Embedder using ProtTrans transformer models from HuggingFace."""

    default_pooling = "mean"

    # Common models for protein embeddings
    PROTEIN_MODELS = {
        "bert": "Rostlab/prot_bert_bfd",
        "prot_bert_bfd": "Rostlab/prot_bert_bfd",
        "t5": "Rostlab/prot_t5_xl_uniref50",
        "albert": "Rostlab/prot_albert",
        "xlnet": "Rostlab/prot_xlnet",
    }

    def __init__(
        self,
        model_name: str = "t5",  # <-- Default is now "t5" instead of "bert"
        layer: int = -1,
        device: Optional[str] = None,
    ):
        """
        Initialize the ProtTrans embedder.

        Args:
            model_name: Model name (t5, bert, albert, xlnet, prot_bert_bfd) or a HuggingFace model path.
                       Default is ProtT5 ("Rostlab/prot_t5_xl_uniref50").
            layer: Which layer to extract embeddings from (-1 for last).
            device: Device to run inference on ('cpu', 'cuda', 'mps', or None for auto).
        """
        # Map short names to full model names
        if model_name in self.PROTEIN_MODELS:
            self.model_name = self.PROTEIN_MODELS[model_name]
            logger.info(f"Using protein model: {self.model_name}")
        else:
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
        logger.info(f"Using device: {self.device} for ProtTrans embeddings")

        # Initialize model
        self._initialize_model(layer)

    def _initialize_model(self, layer):
        """Initialize the transformer model and tokenizer."""
        try:
            from transformers import AutoModel, AutoTokenizer, T5EncoderModel
        except ImportError:
            raise ImportError(
                "Transformers package not found. Install with: pip install transformers"
            )

        logger.info(f"Loading ProtTrans model: {self.model_name}")

        # Load tokenizer
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)

        # Load model based on type
        if "t5" in self.model_name.lower():
            self.model = T5EncoderModel.from_pretrained(self.model_name)
            self.model_type = "t5"
        else:
            self.model = AutoModel.from_pretrained(self.model_name, output_hidden_states=True)
            self.model_type = "bert"  # Default to BERT-like models

        self.model.to(self.device).eval()

        # Count layers
        if self.model_type == "t5":
            self.model_layers = self.model.encoder.config.num_layers
        else:
            self.model_layers = self.model.config.num_hidden_layers

        # Determine which layer to use
        self.repr_layer = self.model_layers + layer if layer < 0 else layer
        if not 0 <= self.repr_layer < self.model_layers:
            raise ValueError(f"Layer {layer} out of bounds (0-{self.model_layers - 1})")

        # Get embedding dimension
        self.embedding_dim = self.model.config.hidden_size

        logger.info(
            f"Model loaded with {sum(p.numel() for p in self.model.parameters()) / 1e6:.1f}M parameters"
        )

    def generate(
        self,
        sequence: str,
        pooling: str = "auto",
        max_length: Optional[int] = 510,  # Account for special tokens
    ) -> np.ndarray:
        """
        Generate transformer embedding for a protein sequence.

        Args:
            sequence: Amino acid sequence
            pooling: Pooling method
            max_length: Maximum sequence length

        Returns:
            Embedding array
        """
        # Handle empty sequence case
        if not sequence:
            if pooling == "none" or (pooling == "auto" and self.default_pooling == "none"):
                return np.zeros((0, self.embedding_dim))
            else:
                return np.zeros(self.embedding_dim)

        # Apply length constraint
        if max_length and len(sequence) > max_length:
            sequence = sequence[:max_length]
            logger.warning(f"Sequence truncated to {max_length} residues")

        # Insert spaces for T5 or BERT
        processed_seq = " ".join(sequence)

        # Tokenize
        inputs = self.tokenizer(
            processed_seq,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=max_length + 2,  # Account for special tokens
        )

        # Generate embeddings
        with torch.no_grad():
            inputs = {k: v.to(self.device) for k, v in inputs.items()}
            outputs = self.model(**inputs)

            if self.model_type == "t5":
                embedding = outputs.last_hidden_state[0].cpu().numpy()
            else:
                embedding = (
                    outputs.hidden_states[self.repr_layer][0, 1 : len(sequence) + 1].cpu().numpy()
                )

        return self._apply_pooling(embedding, pooling)

    def batch_generate(
        self,
        sequences: List[str],
        pooling: str = "auto",
        max_length: Optional[int] = 510,
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

        processed_sequences = []
        seq_lengths = []

        # Prepare sequences
        for seq in sequences:
            if not seq:
                continue

            # Truncate if too long
            if max_length and len(seq) > max_length:
                seq = seq[:max_length]

            processed_sequences.append(" ".join(seq))
            seq_lengths.append(len(seq))

        # Prepare placeholders
        results = []
        for seq in sequences:
            if not seq:
                if pooling == "none" or (pooling == "auto" and self.default_pooling == "none"):
                    results.append(np.zeros((0, self.embedding_dim)))
                else:
                    results.append(np.zeros(self.embedding_dim))
            else:
                results.append(None)

        # Early exit if all empty
        if not processed_sequences:
            return results

        # Tokenize batch
        inputs = self.tokenizer(
            processed_sequences,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=max_length + 2,  # Account for special tokens
        )

        # Generate embeddings
        with torch.no_grad():
            inputs = {k: v.to(self.device) for k, v in inputs.items()}
            outputs = self.model(**inputs)

            if self.model_type == "t5":
                embeddings = outputs.last_hidden_state.cpu().numpy()
                seq_idx = 0
                for i, seq in enumerate(sequences):
                    if not seq:
                        continue
                    emb = embeddings[seq_idx, : seq_lengths[seq_idx]]
                    results[i] = self._apply_pooling(emb, pooling)
                    seq_idx += 1
            else:
                hidden_states = outputs.hidden_states[self.repr_layer].cpu().numpy()
                seq_idx = 0
                for i, seq in enumerate(sequences):
                    if not seq:
                        continue
                    emb = hidden_states[seq_idx, 1 : seq_lengths[seq_idx] + 1]
                    results[i] = self._apply_pooling(emb, pooling)
                    seq_idx += 1

        return results
