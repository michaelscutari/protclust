"""Remote API-based embeddings for protein sequences."""

import hashlib
import json
import os
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import requests

from ..logger import logger
from .baseline import BaseEmbedder


class ESMAPIEmbedder(BaseEmbedder):
    """Embedder using ESM API for protein embeddings."""

    default_pooling = "mean"
    API_URL = "https://api.esmatlas.com/embeddings"

    # ESM model dimensions
    ESM_DIMENSIONS = {
        "esm2_t33_650M_UR50D": 1280,
        "esm2_t36_3B_UR50D": 2560,
        "esm2_t48_15B_UR50D": 5120,
    }

    def __init__(
        self,
        model_name: str = "esm2_t33_650M_UR50D",
        api_key: Optional[str] = None,
        cache_dir: Optional[str] = None,
        max_retries: int = 3,
    ):
        """
        Initialize the ESM API embedder.

        Args:
            model_name: ESM model to use
            api_key: API key for authentication (also checks ESM_API_KEY env var)
            cache_dir: Directory to cache embeddings
            max_retries: Maximum number of retry attempts
        """
        self.model_name = model_name

        # Validate model
        if model_name not in self.ESM_DIMENSIONS:
            raise ValueError(
                f"Unknown ESM model: {model_name}. "
                f"Available models: {list(self.ESM_DIMENSIONS.keys())}"
            )

        # Set embedding dimension
        self.embedding_dim = self.ESM_DIMENSIONS[model_name]

        # Set API key
        self.api_key = api_key or os.environ.get("ESM_API_KEY")
        if not self.api_key:
            logger.warning("No API key provided. Set api_key or ESM_API_KEY environment variable.")

        # Setup caching
        self.use_cache = cache_dir is not None
        if self.use_cache:
            self.cache_dir = Path(cache_dir)
            os.makedirs(self.cache_dir, exist_ok=True)

        self.max_retries = max_retries

    def _get_cache_path(self, sequence: str) -> Path:
        """Get cache file path for a sequence."""
        if not self.use_cache:
            return None

        seq_hash = hashlib.md5(sequence.encode()).hexdigest()
        return self.cache_dir / f"{self.model_name}_{seq_hash}.npz"

    def _get_cached_embedding(self, sequence: str) -> Optional[np.ndarray]:
        """Get embedding from cache if available."""
        if not self.use_cache:
            return None

        cache_path = self._get_cache_path(sequence)
        if cache_path and cache_path.exists():
            try:
                data = np.load(cache_path, allow_pickle=True)
                return data["embedding"]
            except Exception as e:
                logger.warning(f"Failed to load from cache: {e}")
        return None

    def _save_to_cache(self, sequence: str, embedding: np.ndarray) -> None:
        """Save embedding to cache."""
        if not self.use_cache:
            return

        cache_path = self._get_cache_path(sequence)
        try:
            np.savez_compressed(cache_path, embedding=embedding)
        except Exception as e:
            logger.warning(f"Failed to save to cache: {e}")

    def generate(
        self,
        sequence: str,
        pooling: str = "auto",
        max_length: Optional[int] = 1022,
    ) -> np.ndarray:
        """
        Generate ESM embedding via API with caching.

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

        # Try to get from cache
        embedding = self._get_cached_embedding(sequence)
        if embedding is not None:
            logger.debug(f"Using cached embedding for sequence of length {len(sequence)}")
            return self._apply_pooling(embedding, pooling)

        # Prepare API request
        headers = {"Content-Type": "application/json"}
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"

        data = {
            "model": self.model_name,
            "sequence": sequence,
            "representations": ["per_tok"],
        }

        # Make API request with retries
        for attempt in range(self.max_retries):
            try:
                logger.debug(f"Requesting embedding for sequence of length {len(sequence)}")
                response = requests.post(
                    self.API_URL,
                    headers=headers,
                    json=data,
                    timeout=60,
                )
                response.raise_for_status()
                break
            except requests.exceptions.RequestException as e:
                logger.warning(
                    f"API request failed (attempt {attempt + 1}/{self.max_retries}): {e}"
                )
                if attempt < self.max_retries - 1:
                    time.sleep(2**attempt)  # Exponential backoff
                else:
                    raise ValueError(
                        f"Failed to get embedding from API after {self.max_retries} attempts"
                    ) from e

        # Process response
        try:
            result = response.json()

            # Extract embeddings (format depends on API)
            if "embeddings" in result:
                embedding = np.array(result["embeddings"])
            elif "representations" in result and "per_tok" in result["representations"]:
                embedding = np.array(result["representations"]["per_tok"])
            else:
                raise ValueError("Unexpected API response format")

            # Cache result
            self._save_to_cache(sequence, embedding)

            # Apply pooling
            return self._apply_pooling(embedding, pooling)

        except (json.JSONDecodeError, KeyError) as e:
            raise ValueError(f"Failed to parse API response: {e}")

    def batch_generate(
        self,
        sequences: List[str],
        pooling: str = "auto",
        max_length: Optional[int] = 1022,
    ) -> List[np.ndarray]:
        """
        Generate embeddings for multiple sequences.

        Args:
            sequences: List of amino acid sequences
            pooling: Pooling method
            max_length: Maximum sequence length

        Returns:
            List of embedding arrays
        """
        # API doesn't support true batching, so we process sequentially
        results = []

        for sequence in sequences:
            results.append(self.generate(sequence, pooling, max_length))

        return results
