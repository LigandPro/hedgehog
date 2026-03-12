"""Tests for remote Matcha checkout selection."""

import pytest

import hedgehog.setup as setup_pkg
from hedgehog.setup._matcha import _select_matcha_pr


def _patch_matcha_api(monkeypatch, pulls_payload, reviews_payload):
    def _fake_query(url: str):
        if "pulls?state=open" in url:
            return pulls_payload
        for pr_number, payload in reviews_payload.items():
            if f"/pulls/{pr_number}/reviews" in url:
                return payload
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("hedgehog.setup._matcha._query_json", _fake_query)


def test_select_matcha_pr_prefers_latest_unapproved(monkeypatch):
    pulls_payload = [
        {
            "number": 40,
            "draft": False,
            "title": "latest",
            "html_url": "https://example.test/40",
            "updated_at": "2026-03-11T07:31:30Z",
            "head": {"sha": "sha40", "ref": "pr40"},
        },
        {
            "number": 39,
            "draft": False,
            "title": "older",
            "html_url": "https://example.test/39",
            "updated_at": "2026-03-10T07:31:30Z",
            "head": {"sha": "sha39", "ref": "pr39"},
        },
    ]
    reviews_payload = {
        40: [],
        39: [{"state": "COMMENTED"}],
    }

    _patch_matcha_api(monkeypatch, pulls_payload, reviews_payload)

    selected = _select_matcha_pr()
    assert selected["number"] == 40
    assert selected["head"]["sha"] == "sha40"


def test_select_matcha_pr_skips_approved_pr(monkeypatch):
    pulls_payload = [
        {
            "number": 40,
            "draft": False,
            "title": "approved",
            "html_url": "https://example.test/40",
            "updated_at": "2026-03-11T07:31:30Z",
            "head": {"sha": "sha40", "ref": "pr40"},
        },
        {
            "number": 39,
            "draft": False,
            "title": "unapproved",
            "html_url": "https://example.test/39",
            "updated_at": "2026-03-10T07:31:30Z",
            "head": {"sha": "sha39", "ref": "pr39"},
        },
    ]

    _patch_matcha_api(
        monkeypatch,
        pulls_payload,
        {
            40: [{"state": "APPROVED"}],
            39: [],
        },
    )

    selected = _select_matcha_pr()
    assert selected["number"] == 39
    assert selected["head"]["sha"] == "sha39"


def test_select_matcha_pr_prefers_higher_number_on_equal_updated_at(monkeypatch):
    pulls_payload = [
        {
            "number": 38,
            "draft": False,
            "title": "older-number",
            "html_url": "https://example.test/38",
            "updated_at": "2026-03-11T07:31:30Z",
            "head": {"sha": "sha38", "ref": "pr38"},
        },
        {
            "number": 40,
            "draft": False,
            "title": "newer-number",
            "html_url": "https://example.test/40",
            "updated_at": "2026-03-11T07:31:30Z",
            "head": {"sha": "sha40", "ref": "pr40"},
        },
    ]

    _patch_matcha_api(monkeypatch, pulls_payload, {38: [], 40: []})

    selected = _select_matcha_pr()
    assert selected["number"] == 40
    assert selected["head"]["sha"] == "sha40"


def test_select_matcha_pr_raises_when_all_open_prs_are_approved(monkeypatch):
    pulls_payload = [
        {
            "number": 40,
            "draft": False,
            "title": "approved-a",
            "html_url": "https://example.test/40",
            "updated_at": "2026-03-11T07:31:30Z",
            "head": {"sha": "sha40", "ref": "pr40"},
        },
        {
            "number": 39,
            "draft": False,
            "title": "approved-b",
            "html_url": "https://example.test/39",
            "updated_at": "2026-03-10T07:31:30Z",
            "head": {"sha": "sha39", "ref": "pr39"},
        },
    ]

    _patch_matcha_api(
        monkeypatch,
        pulls_payload,
        {
            40: [{"state": "APPROVED"}],
            39: [{"state": "APPROVED"}],
        },
    )

    with pytest.raises(RuntimeError, match="No eligible Matcha pull request was found"):
        _select_matcha_pr()


def test_setup_dir_exposes_other_lazy_exports_without_binding():
    setup_pkg.__dict__.pop("ensure_nvmolkit_worker", None)

    assert "ensure_nvmolkit_worker" in dir(setup_pkg)
    assert "ensure_nvmolkit_worker" not in setup_pkg.__dict__
    assert callable(setup_pkg.ensure_nvmolkit_worker)
