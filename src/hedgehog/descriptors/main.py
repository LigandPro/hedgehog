"""Compatibility wrapper for the descriptors stage."""

from hedgehog.descriptors.stage import run


def main(data, config, subfolder=None, reporter=None):
    """Run the canonical descriptors stage entry point."""
    return run(data, config, subfolder=subfolder, reporter=reporter)
