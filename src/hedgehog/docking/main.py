"""Backward-compatible docking entry point."""

from __future__ import annotations

from hedgehog.docking.stage import run


def main(config: dict, reporter=None):
    """Run docking stage via the current stage entry point."""
    return run(config, reporter=reporter)
