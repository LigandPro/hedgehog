"""Focused tests for shared path helpers and lazy package exports."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

from hedgehog.descriptors.io import process_path as descriptors_process_path
from hedgehog.struct_filters.io import process_path as struct_filters_io_process_path
from hedgehog.utils.paths import process_path


def _reload_package_without_submodules(
    monkeypatch, package_name: str, submodules: list[str]
):
    package = importlib.import_module(package_name)
    for module_name in submodules:
        parent_name, _, child_name = module_name.rpartition(".")
        parent_module = importlib.import_module(parent_name)
        monkeypatch.delattr(parent_module, child_name, raising=False)
        monkeypatch.delitem(sys.modules, module_name, raising=False)
    return importlib.reload(package)


def test_shared_process_path_accepts_path_objects(tmp_path):
    result = process_path(tmp_path / "out")
    assert result.endswith("/")
    assert Path(result.rstrip("/")).exists()


def test_shared_process_path_appends_keyword_subfolder(tmp_path):
    result = process_path(tmp_path, "subdir")
    assert result.endswith("subdir/")
    assert Path(result.rstrip("/")).exists()


def test_struct_filters_io_process_path_stays_compatible(tmp_path):
    target = tmp_path / "compat_out"
    assert struct_filters_io_process_path(target) == process_path(target)


def test_descriptors_io_reexports_shared_process_path():
    assert descriptors_process_path is process_path


def test_docking_package_exports_run_lazily(monkeypatch):
    module = _reload_package_without_submodules(
        monkeypatch,
        "hedgehog.docking",
        ["hedgehog.docking.stage"],
    )
    assert "hedgehog.docking.stage" not in sys.modules
    assert callable(module.run)
    assert "hedgehog.docking.stage" in sys.modules


def test_descriptors_package_exports_run_lazily(monkeypatch):
    module = _reload_package_without_submodules(
        monkeypatch,
        "hedgehog.descriptors",
        ["hedgehog.descriptors.stage"],
    )
    assert "hedgehog.descriptors.stage" not in sys.modules
    assert callable(module.run)
    assert "hedgehog.descriptors.stage" in sys.modules


def test_reporting_package_exports_report_generator_lazily(monkeypatch):
    module = _reload_package_without_submodules(
        monkeypatch,
        "hedgehog.reporting",
        ["hedgehog.reporting.report_generator"],
    )
    assert "hedgehog.reporting.report_generator" not in sys.modules
    assert module.ReportGenerator.__name__ == "ReportGenerator"
    assert "hedgehog.reporting.report_generator" in sys.modules


def test_docking_filters_package_exports_main_lazily(monkeypatch):
    module = _reload_package_without_submodules(
        monkeypatch,
        "hedgehog.docking_filters",
        [
            "hedgehog.docking_filters.main",
            "hedgehog.docking_filters.utils",
        ],
    )
    assert "hedgehog.docking_filters.main" not in sys.modules
    assert "hedgehog.docking_filters.utils" not in sys.modules
    assert callable(module.docking_filters_main)
    assert "hedgehog.docking_filters.main" in sys.modules
