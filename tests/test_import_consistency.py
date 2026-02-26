"""Static import consistency checks for internal modules."""

from __future__ import annotations

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"

# Optional integration module used behind a guarded import.
OPTIONAL_MISSING_MODULES = {
    "hedgehog.vendor.moleval.metrics.posecheck",
}


def _iter_python_files() -> list[Path]:
    files: list[Path] = []
    for root in (PROJECT_ROOT / "src" / "hedgehog", PROJECT_ROOT / "tests"):
        files.extend(root.rglob("*.py"))
    return files


def _collect_hedgehog_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "hedgehog" or alias.name.startswith("hedgehog."):
                    modules.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level != 0 or not node.module:
                continue
            if node.module == "hedgehog" or node.module.startswith("hedgehog."):
                modules.add(node.module)
    return modules


def _module_exists(module_name: str) -> bool:
    module_path = SRC_ROOT / module_name.replace(".", "/")
    return module_path.with_suffix(".py").exists() or module_path.exists()


def test_internal_hedgehog_imports_point_to_existing_modules() -> None:
    missing: dict[str, list[str]] = {}
    for path in _iter_python_files():
        rel_path = path.relative_to(PROJECT_ROOT)
        for module_name in _collect_hedgehog_imports(path):
            if module_name in OPTIONAL_MISSING_MODULES:
                continue
            if not _module_exists(module_name):
                missing.setdefault(module_name, []).append(str(rel_path))

    assert not missing, "Found imports to missing internal modules:\n" + "\n".join(
        f"- {module_name} (imported in {', '.join(sorted(files))})"
        for module_name, files in sorted(missing.items())
    )
