"""Compatibility tests for docking entry points."""

from __future__ import annotations


def test_docking_main_delegates_to_stage_run(monkeypatch) -> None:
    import hedgehog.docking.main as docking_main

    called: dict[str, object] = {}

    def _fake_run(config, reporter=None):
        called["config"] = config
        called["reporter"] = reporter
        return "ok"

    monkeypatch.setattr(docking_main, "run", _fake_run)
    result = docking_main.main({"k": "v"}, reporter="r")

    assert result == "ok"
    assert called == {"config": {"k": "v"}, "reporter": "r"}
