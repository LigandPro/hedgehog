"""Auto-install/update Matcha checkout from the latest open unapproved PR."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Any

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download

_MATCHA_REPO_URL = "https://github.com/LigandPro/Matcha.git"
_MATCHA_API_PULLS = (
    "https://api.github.com/repos/LigandPro/Matcha/pulls"
    "?state=open&sort=updated&direction=desc&per_page=20"
)


def _github_request(url: str) -> urllib.request.Request:
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "hedgehog-matcha-bootstrap",
    }
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return urllib.request.Request(url, headers=headers)


def _query_json(url: str) -> Any:
    with urllib.request.urlopen(_github_request(url), timeout=30) as response:  # noqa: S310
        return json.loads(response.read())


def _pr_reviews_url(pr_number: int) -> str:
    return (
        f"https://api.github.com/repos/LigandPro/Matcha/pulls/{pr_number}/reviews"
        "?per_page=100"
    )


def _pr_has_approval(pr_number: int) -> bool:
    reviews = _query_json(_pr_reviews_url(pr_number))
    return any(str(review.get("state", "")).upper() == "APPROVED" for review in reviews)


def _select_matcha_pr(pr_number: int | None = None) -> dict[str, Any]:
    pulls = _query_json(_MATCHA_API_PULLS)
    if not isinstance(pulls, list) or not pulls:
        raise RuntimeError("No open Matcha pull requests were returned by GitHub.")

    pulls = sorted(
        pulls,
        key=lambda pr: (str(pr.get("updated_at", "")), int(pr.get("number", 0))),
        reverse=True,
    )
    candidates = [pr for pr in pulls if not pr.get("draft")]
    fallback: dict[str, Any] | None = candidates[0] if candidates else None

    review_lookup_failed = False
    for pr in candidates:
        number = int(pr["number"])
        if pr_number is not None and number != pr_number:
            continue
        try:
            has_approval = _pr_has_approval(number)
        except Exception as exc:  # noqa: BLE001
            review_lookup_failed = True
            logger.warning(
                "Could not load Matcha PR #%d reviews from GitHub: %s", number, exc
            )
            continue
        if not has_approval:
            return pr

    if pr_number is not None:
        raise RuntimeError(
            f"Matcha PR #{pr_number} is not open or already has an approval."
        )
    if review_lookup_failed and fallback is not None:
        logger.warning(
            "Falling back to latest open Matcha PR #%s%s.",
            fallback["number"],
            " because review lookup failed." if review_lookup_failed else "",
        )
        return fallback
    raise RuntimeError("No eligible Matcha pull request was found.")


def _run_git(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, timeout=1800)


def _ensure_checkout_parent(checkout_dir: Path) -> None:
    checkout_dir.parent.mkdir(parents=True, exist_ok=True)


def _clone_or_refresh_repo(checkout_dir: Path, repo_url: str) -> None:
    if checkout_dir.exists() and not (checkout_dir / ".git").exists():
        logger.warning(
            "Removing incomplete Matcha checkout at %s before recloning.", checkout_dir
        )
        shutil.rmtree(checkout_dir)

    if not checkout_dir.exists():
        _ensure_checkout_parent(checkout_dir)
        _run_git(["git", "clone", repo_url, str(checkout_dir)], cwd=checkout_dir.parent)
        return

    _run_git(["git", "fetch", "--all", "--prune"], cwd=checkout_dir)


def _checkout_pr_head(checkout_dir: Path, pr: dict[str, Any]) -> None:
    pr_number = int(pr["number"])
    local_branch = f"hedgehog-matcha-pr-{pr_number}"
    _run_git(
        ["git", "fetch", "--force", "origin", f"pull/{pr_number}/head:{local_branch}"],
        cwd=checkout_dir,
    )
    _run_git(["git", "checkout", "--detach", str(pr["head"]["sha"])], cwd=checkout_dir)
    metadata = {
        "pr_number": pr_number,
        "title": pr.get("title"),
        "url": pr.get("html_url"),
        "head_sha": pr["head"]["sha"],
        "head_ref": pr["head"]["ref"],
        "updated_at": pr.get("updated_at"),
    }
    (checkout_dir / ".hedgehog_matcha_pr.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def _has_matcha_cli(checkout_dir: Path) -> bool:
    """Return True when the expected Matcha CLI module exists."""
    return (checkout_dir / "matcha" / "cli.py").exists()


def _checkout_main_head(checkout_dir: Path) -> None:
    """Fallback to repository main branch HEAD."""
    _run_git(["git", "fetch", "origin", "main"], cwd=checkout_dir)
    _run_git(["git", "checkout", "--detach", "origin/main"], cwd=checkout_dir)
    metadata = {
        "source": "origin/main",
        "note": "fallback because selected PR checkout lacked matcha/cli.py",
    }
    (checkout_dir / ".hedgehog_matcha_pr.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def ensure_matcha_checkout(
    project_root: Path,
    *,
    repo_url: str = _MATCHA_REPO_URL,
    checkout_dir: str | Path | None = None,
    pr_number: int | None = None,
) -> Path:
    """Ensure a managed Matcha checkout exists and return its path."""
    if not shutil.which("git"):
        raise RuntimeError("git is required to fetch Matcha from GitHub.")

    target_dir = (
        Path(checkout_dir).expanduser()
        if checkout_dir is not None
        else (project_root / "modules" / "matcha_remote")
    )
    if not target_dir.is_absolute():
        target_dir = (project_root / target_dir).resolve()

    pr = _select_matcha_pr(pr_number=pr_number)
    size_hint = "~repository checkout plus uv environment for Matcha"
    if not target_dir.exists() and not confirm_download(
        "Matcha remote checkout", size_hint
    ):
        raise RuntimeError("Matcha download declined by user.")

    logger.info(
        "Preparing Matcha checkout from PR #%s (%s)",
        pr["number"],
        pr["head"]["sha"],
    )
    _clone_or_refresh_repo(target_dir, repo_url)
    _checkout_pr_head(target_dir, pr)
    if not _has_matcha_cli(target_dir):
        logger.warning(
            "Selected Matcha PR checkout is missing matcha/cli.py. Falling back to origin/main."
        )
        _checkout_main_head(target_dir)
        if not _has_matcha_cli(target_dir):
            raise RuntimeError(
                "Matcha checkout is missing matcha/cli.py after fallback to origin/main."
            )
    return target_dir
