"""Backward-compatible docking stage entry point."""

from __future__ import annotations

from hedgehog.docking.main import main


def run(config: dict, reporter=None):
    """Run docking stage via the current main entry point."""
    return main(config, reporter=reporter)
