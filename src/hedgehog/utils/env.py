"""Environment variable helpers."""

import os


def plain_output_enabled() -> bool:
    """Check if plain (non-Rich) output is requested via HEDGEHOG_PLAIN_OUTPUT=1."""
    return os.environ.get("HEDGEHOG_PLAIN_OUTPUT", "").strip() == "1"
