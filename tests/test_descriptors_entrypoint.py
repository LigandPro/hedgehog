"""Tests for descriptors stage entry points."""

import pandas as pd


def test_descriptors_main_delegates_to_stage_run(monkeypatch):
    """Legacy descriptors.main should use the canonical stage implementation."""
    import hedgehog.descriptors.main as descriptors_main

    captured = {}

    def _fake_run(data, config, subfolder=None, reporter=None):
        captured["data"] = data
        captured["config"] = config
        captured["subfolder"] = subfolder
        captured["reporter"] = reporter
        return "ok"

    monkeypatch.setattr(descriptors_main, "run", _fake_run)

    data = pd.DataFrame({"smiles": ["CCO"]})
    config = {"config_descriptors": "config.yml"}
    reporter = object()

    assert (
        descriptors_main.main(data, config, subfolder="custom", reporter=reporter)
        == "ok"
    )
    assert captured["data"] is data
    assert captured["config"] is config
    assert captured["subfolder"] == "custom"
    assert captured["reporter"] is reporter
