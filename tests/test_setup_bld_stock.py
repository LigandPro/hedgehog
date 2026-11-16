"""Tests for the optional BLDpharm AiZynthFinder stock (no network access)."""

from __future__ import annotations

import csv
import io
import json
import urllib.parse
from pathlib import Path

import pytest
import yaml
from rdkit import Chem

from hedgehog.setup import _bld_stock
from hedgehog.setup._bld_stock import (
    build_bld_stock,
    ensure_bld_stock,
    get_bld_stock_path,
    register_stock,
)

# Tiny stand-in for the PubChem deposit: CID 2 and 3 share one structure,
# CID 4 cannot be parsed, CID 1 is listed twice by two substance records.
_SMILES = {
    1: "OB(O)c1ccccc1",
    2: "CCO",
    3: "OCC",
    4: "C1CC",
    5: "c1ccncc1",
}


class FakePubChem:
    """Serve the fixture through the PubChem request interface."""

    def __init__(self) -> None:
        self.smiles_requests: list[list[int]] = []

    def __call__(self, url: str, data: bytes | None) -> bytes:
        if data is None:
            assert url.endswith("/substance/sourceall/BLD%20Pharm/cids/TXT")
            return b"5\n1\n2\n3\n4\n1\n"
        assert url.endswith("/compound/cid/property/SMILES/CSV")
        cids = [
            int(c) for c in urllib.parse.parse_qs(data.decode())["cid"][0].split(",")
        ]
        self.smiles_requests.append(cids)
        out = io.StringIO()
        writer = csv.writer(out, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(["CID", "SMILES"])
        writer.writerows([cid, _SMILES[cid]] for cid in cids)
        return out.getvalue().encode()


def _key(smiles: str) -> str:
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


def _stock_keys(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "inchi_key"
    return lines[1:]


def test_build_writes_rdkit_inchikey_stock_and_manifest(tmp_path: Path):
    stock = tmp_path / "bld" / "bld_stock.csv"
    stock.parent.mkdir()
    fake = FakePubChem()

    manifest = build_bld_stock(stock, fetch=fake, batch_size=2)

    expected = sorted({_key(_SMILES[c]) for c in (1, 2, 5)})
    assert _stock_keys(stock) == expected
    assert sorted(fake.smiles_requests) == [[1, 2], [3, 4], [5]]
    assert manifest["pubchem_source_id"] == 22134
    assert manifest["unique_cids"] == 5
    assert manifest["invalid_cids"] == [4]
    assert manifest["molecules"] == 3
    assert "not provided" in manifest["availability"]
    assert json.loads((stock.parent / "manifest.json").read_text()) == manifest
    assert not (stock.parent / ".build").exists()


def test_build_resumes_from_cached_pages(tmp_path: Path):
    stock = tmp_path / "bld_stock.csv"
    pages = tmp_path / ".build" / "pages"
    pages.mkdir(parents=True)
    (tmp_path / ".build" / "cids.txt").write_text("1\n2\n3\n")
    (pages / "00000.csv").write_text('"CID","SMILES"\n1,"OB(O)c1ccccc1"\n2,"CCO"\n')
    fake = FakePubChem()

    build_bld_stock(stock, fetch=fake, batch_size=2)

    assert fake.smiles_requests == [[3]]
    assert _stock_keys(stock) == sorted({_key("OB(O)c1ccccc1"), _key("CCO")})


def test_build_rejects_incomplete_pubchem_batch(tmp_path: Path):
    stock = tmp_path / "bld_stock.csv"

    def truncated(url: str, data: bytes | None) -> bytes:
        if data is None:
            return b"1\n2\n"
        return b'"CID","SMILES"\n1,"CCO"\n'

    with pytest.raises(RuntimeError, match="returned 1 of 2"):
        build_bld_stock(stock, fetch=truncated)
    assert not stock.exists()


def test_register_stock_keeps_existing_entries(tmp_path: Path):
    config = tmp_path / "config.yml"
    config.write_text(
        yaml.safe_dump(
            {
                "expansion": {"uspto": ["m.onnx", "t.csv.gz"]},
                "stock": {"zinc": "zinc.hdf5"},
            }
        )
    )

    register_stock(config, "bld", tmp_path / "bld_stock.csv")
    register_stock(config, "bld", tmp_path / "bld_stock.csv")

    data = yaml.safe_load(config.read_text())
    assert data["stock"] == {
        "zinc": "zinc.hdf5",
        "bld": str(tmp_path / "bld_stock.csv"),
    }
    assert data["expansion"] == {"uspto": ["m.onnx", "t.csv.gz"]}


def test_ensure_reuses_existing_stock_and_registers_it(tmp_path: Path, monkeypatch):
    config = tmp_path / "modules" / "aizynthfinder" / "public" / "config.yml"
    config.parent.mkdir(parents=True)
    config.write_text("stock:\n  zinc: zinc.hdf5\n")
    stock = get_bld_stock_path(tmp_path)
    stock.parent.mkdir(parents=True)
    stock.write_text("inchi_key\nLFQSCWFLJHTTHZ-UHFFFAOYSA-N\n")
    (stock.parent / "manifest.json").write_text("{}\n")
    monkeypatch.setattr(
        "hedgehog.setup._aizynthfinder.ensure_aizynthfinder", lambda root: config
    )

    def unexpected_build(*args, **kwargs):
        raise AssertionError("existing stock must not be rebuilt")

    monkeypatch.setattr(_bld_stock, "build_bld_stock", unexpected_build)

    assert ensure_bld_stock(tmp_path) == stock
    assert yaml.safe_load(config.read_text())["stock"]["bld"] == str(stock)


def test_ensure_respects_declined_download(tmp_path: Path, monkeypatch):
    config = tmp_path / "config.yml"
    config.write_text("stock: {}\n")
    monkeypatch.setattr(
        "hedgehog.setup._aizynthfinder.ensure_aizynthfinder", lambda root: config
    )
    monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
    monkeypatch.setenv("HEDGEHOG_NON_INTERACTIVE", "1")

    with pytest.raises(RuntimeError, match="declined"):
        ensure_bld_stock(tmp_path)
    assert "bld" not in (yaml.safe_load(config.read_text())["stock"] or {})
