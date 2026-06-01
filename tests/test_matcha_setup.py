"""Tests for Matcha checkout bootstrap."""

import hedgehog.setup as setup_pkg
from hedgehog.setup._matcha import ensure_matcha_checkout


def test_ensure_matcha_checkout_clones_and_checks_out_main(monkeypatch, tmp_path):
    """Clone LigandPro/Matcha and detach at origin/main."""
    target = tmp_path / "modules" / "matcha_remote"
    calls: list[list[str]] = []

    def run_git(cmd, cwd):
        calls.append(list(cmd))
        if cmd[:2] == ["git", "clone"]:
            target.mkdir(parents=True, exist_ok=True)
            (target / ".git").mkdir()
            (target / "matcha").mkdir()
            (target / "matcha" / "cli.py").write_text("# stub\n", encoding="utf-8")

    monkeypatch.setattr("hedgehog.setup._matcha._run_git", run_git)
    monkeypatch.setattr("hedgehog.setup._matcha.shutil.which", lambda _: "/usr/bin/git")
    monkeypatch.setattr("hedgehog.setup._matcha.confirm_download", lambda *a: True)

    result = ensure_matcha_checkout(tmp_path)

    assert result == target.resolve()
    assert [
        "git",
        "clone",
        "https://github.com/LigandPro/Matcha.git",
        str(target),
    ] in calls
    assert ["git", "checkout", "--detach", "origin/main"] in calls
    assert not any(cmd[:3] == ["git", "fetch", "origin"] for cmd in calls)


def test_ensure_matcha_checkout_updates_existing_repo(monkeypatch, tmp_path):
    """An existing checkout is fetched and reset to origin/main."""
    target = tmp_path / "modules" / "matcha_remote"
    target.mkdir(parents=True)
    (target / ".git").mkdir()
    (target / "matcha").mkdir()
    (target / "matcha" / "cli.py").write_text("# stub\n", encoding="utf-8")

    calls: list[list[str]] = []

    monkeypatch.setattr(
        "hedgehog.setup._matcha._run_git",
        lambda cmd, cwd: calls.append(list(cmd)),
    )
    monkeypatch.setattr("hedgehog.setup._matcha.shutil.which", lambda _: "/usr/bin/git")
    monkeypatch.setattr("hedgehog.setup._matcha.confirm_download", lambda *a: True)

    result = ensure_matcha_checkout(tmp_path)

    assert result == target.resolve()
    assert ["git", "fetch", "origin", "main"] in calls
    assert ["git", "checkout", "--detach", "origin/main"] in calls
    assert not any(cmd[:2] == ["git", "clone"] for cmd in calls)


def test_setup_dir_exposes_other_lazy_exports_without_binding():
    setup_pkg.__dict__.pop("ensure_nvmolkit_worker", None)

    assert "ensure_nvmolkit_worker" in dir(setup_pkg)
    assert "ensure_nvmolkit_worker" not in setup_pkg.__dict__
    assert callable(setup_pkg.ensure_nvmolkit_worker)
