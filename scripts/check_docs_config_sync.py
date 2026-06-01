"""Check documentation for known stale configuration defaults."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STALE_PATTERNS = {
    "n_jobs: 62": "Use n_jobs: -1 from src/hedgehog/configs/config.yml.",
    "n_jobs: 4": "Use n_jobs: -1 from config_moleval.yml unless documenting an example.",
    "sa_score_min: 0": "Use sa_score_min: 1 from config_synthesis.yml.",
    "mce18_min: 45": "Use mce18_min: 20 from config_descriptors.yml.",
    "mce18_max: 100": "Use mce18_max: 140 from config_descriptors.yml.",
    "nibr_scheduler: threads": "Use nibr_scheduler: processes from config_structFilters.yml.",
    "required_residues: []": "Use required_residues: ['ASP12'] from config_docking_filters.yml.",
}


def main() -> int:
    failures: list[str] = []
    for path in sorted((ROOT / "docs" / "content").rglob("*.mdx")):
        text = path.read_text(encoding="utf-8")
        rel = path.relative_to(ROOT)
        for pattern, message in STALE_PATTERNS.items():
            if pattern in text:
                failures.append(f"{rel}: contains {pattern!r}. {message}")

    if failures:
        print("Documentation/config sync check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print("Documentation/config sync check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
