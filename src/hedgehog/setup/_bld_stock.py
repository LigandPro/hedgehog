"""Optional BLDpharm building-block stock for AiZynthFinder.

The stock is built from the public PubChem deposit of the supplier BLD Pharm
(PubChem source 22134). It records catalogue membership only: PubChem does not
provide current inventory, prices, or delivery times.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import shutil
import time
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import yaml

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download

BLD_STOCK_KEY = "bld"
BLD_SOURCE_NAME = "BLD Pharm"
BLD_SOURCE_ID = 22134
BLD_SOURCE_URL = f"https://pubchem.ncbi.nlm.nih.gov/source/{BLD_SOURCE_ID}"
BLD_AVAILABILITY_NOTE = (
    "Catalogue membership from PubChem only; current stock, price and "
    "delivery are not provided"
)

_PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_CIDS_URL = (
    f"{_PUBCHEM_REST}/substance/sourceall/"
    f"{urllib.parse.quote(BLD_SOURCE_NAME)}/cids/TXT"
)
_SMILES_URL = f"{_PUBCHEM_REST}/compound/cid/property/SMILES/CSV"
_BATCH_SIZE = 5000
_DOWNLOAD_WORKERS = 3
_RETRIES = 5
_REQUEST_TIMEOUT = 300

Fetch = Callable[[str, bytes | None], bytes]


def get_bld_stock_dir(project_root: Path) -> Path:
    """Return the project-local directory holding the BLDpharm stock."""
    return project_root / "modules" / "aizynthfinder" / BLD_STOCK_KEY


def get_bld_stock_path(project_root: Path) -> Path:
    """Return the AiZynthFinder stock file path for BLDpharm."""
    return get_bld_stock_dir(project_root) / "bld_stock.csv"


def _pubchem_fetch(url: str, data: bytes | None = None) -> bytes:
    """Fetch a PubChem PUG REST response with retries."""
    for attempt in range(_RETRIES):
        try:
            request = urllib.request.Request(url, data=data)  # noqa: S310
            with urllib.request.urlopen(  # noqa: S310
                request, timeout=_REQUEST_TIMEOUT
            ) as response:
                return response.read()
        except OSError as exc:
            if attempt == _RETRIES - 1:
                raise
            delay = 5 * (attempt + 1)
            logger.warning("PubChem request failed (%s); retrying in %ds", exc, delay)
            time.sleep(delay)
    raise AssertionError("unreachable")


def _read_cids(raw: bytes) -> list[int]:
    """Parse PubChem TXT output into a sorted list of unique CIDs."""
    return sorted({int(line) for line in raw.decode().split() if line.strip()})


def _fetch_page(cids: list[int], dest: Path, fetch: Fetch) -> None:
    """Download SMILES for one CID batch into ``dest`` (skipped if present)."""
    if dest.exists():
        return
    body = urllib.parse.urlencode({"cid": ",".join(map(str, cids))}).encode()
    raw = fetch(_SMILES_URL, body)
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    returned = {int(row["CID"]) for row in rows}
    if returned != set(cids):
        raise RuntimeError(
            f"PubChem returned {len(returned)} of {len(cids)} requested CIDs "
            f"for batch starting at CID {cids[0]}"
        )
    tmp = dest.with_suffix(".tmp")
    tmp.write_bytes(raw)
    tmp.replace(dest)


def _page_inchikeys(path: Path) -> tuple[set[str], list[int]]:
    """Return RDKit InChIKeys and unparseable CIDs for one downloaded page."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    keys: set[str] = set()
    invalid: list[int] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            mol = Chem.MolFromSmiles(row["SMILES"] or "")
            key = Chem.MolToInchiKey(mol) if mol is not None else ""
            if key:
                keys.add(key)
            else:
                invalid.append(int(row["CID"]))
    return keys, invalid


def _collect_inchikeys(pages: list[Path], jobs: int) -> tuple[set[str], list[int]]:
    """Convert downloaded pages to InChIKeys, optionally in parallel."""
    keys: set[str] = set()
    invalid: list[int] = []
    results: Iterable[tuple[set[str], list[int]]]
    if jobs > 1 and len(pages) > 1:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(_page_inchikeys, pages))
    else:
        results = map(_page_inchikeys, pages)
    for page_keys, page_invalid in results:
        keys.update(page_keys)
        invalid.extend(page_invalid)
    return keys, sorted(invalid)


def build_bld_stock(
    stock_path: Path,
    *,
    fetch: Fetch = _pubchem_fetch,
    batch_size: int = _BATCH_SIZE,
    jobs: int = 1,
) -> dict:
    """Build the BLDpharm InChIKey stock from PubChem.

    Downloads the CIDs deposited by BLD Pharm, fetches their SMILES in batches,
    and derives InChIKeys with RDKit so that they match the keys AiZynthFinder
    computes for precursors. Partial downloads are cached in ``.build/`` next to
    the stock file, so an interrupted build resumes where it stopped.

    Args:
        stock_path: Destination CSV with an ``inchi_key`` column.
        fetch: Callable performing a PubChem request (injected in tests).
        batch_size: Number of CIDs per SMILES request.
        jobs: Worker processes used for the RDKit conversion.

    Returns:
        The manifest written next to the stock file.
    """
    work_dir = stock_path.parent / ".build"
    pages_dir = work_dir / "pages"
    pages_dir.mkdir(parents=True, exist_ok=True)

    cids_file = work_dir / "cids.txt"
    if not cids_file.exists():
        logger.info("Downloading BLD Pharm CID list from PubChem...")
        cids_file.write_text(
            "\n".join(map(str, _read_cids(fetch(_CIDS_URL, None)))) + "\n",
            encoding="utf-8",
        )
    cids = _read_cids(cids_file.read_bytes())
    if not cids:
        raise RuntimeError("PubChem returned no CIDs for BLD Pharm")

    batches = [cids[i : i + batch_size] for i in range(0, len(cids), batch_size)]
    pages = [pages_dir / f"{index:05d}.csv" for index in range(len(batches))]
    logger.info(
        "Fetching structures for %d BLD Pharm CIDs in %d batches...",
        len(cids),
        len(batches),
    )
    with ThreadPoolExecutor(max_workers=_DOWNLOAD_WORKERS) as pool:
        list(
            pool.map(
                lambda item: _fetch_page(item[0], item[1], fetch),
                zip(batches, pages, strict=True),
            )
        )

    logger.info("Computing RDKit InChIKeys...")
    keys, invalid = _collect_inchikeys(pages, jobs)
    if not keys:
        raise RuntimeError("No valid structures in the BLD Pharm download")

    tmp = stock_path.with_suffix(".csv.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write("inchi_key\n")
        handle.writelines(f"{key}\n" for key in sorted(keys))
    tmp.replace(stock_path)

    manifest = {
        "stock_key": BLD_STOCK_KEY,
        "supplier": BLD_SOURCE_NAME,
        "pubchem_source_id": BLD_SOURCE_ID,
        "source_url": BLD_SOURCE_URL,
        "cid_list_url": _CIDS_URL,
        "retrieved_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "availability": BLD_AVAILABILITY_NOTE,
        "unique_cids": len(cids),
        "invalid_cids": invalid,
        "molecules": len(keys),
        "inchikey_method": "RDKit MolToInchiKey from PubChem SMILES",
        "cid_list_sha256": hashlib.sha256(cids_file.read_bytes()).hexdigest(),
        "stock_sha256": hashlib.sha256(stock_path.read_bytes()).hexdigest(),
    }
    _manifest_path(stock_path).write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    shutil.rmtree(work_dir)
    logger.info(
        "BLDpharm stock: %d unique InChIKeys from %d CIDs (%d unparseable)",
        len(keys),
        len(cids),
        len(invalid),
    )
    return manifest


def _manifest_path(stock_path: Path) -> Path:
    return stock_path.with_name("manifest.json")


def register_stock(config_yml: Path, key: str, stock_path: Path) -> None:
    """Add or update one stock entry in an AiZynthFinder ``config.yml``."""
    with config_yml.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    stock = config.get("stock") or {}
    if stock.get(key) == str(stock_path):
        return
    stock[key] = str(stock_path)
    config["stock"] = stock
    tmp = config_yml.with_suffix(".yml.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    tmp.replace(config_yml)
    logger.info("Registered AiZynthFinder stock '%s' in %s", key, config_yml)


def ensure_bld_stock(
    project_root: Path, *, refresh: bool = False, jobs: int = 1
) -> Path:
    """Ensure the BLDpharm stock exists and is registered with AiZynthFinder.

    Installs AiZynthFinder first when needed, builds
    ``modules/aizynthfinder/bld/bld_stock.csv`` from PubChem unless it already
    exists (or ``refresh`` is set), and adds it to ``config.yml`` under the
    ``bld`` key. Once registered, AiZynthFinder searches it together with the
    other configured stocks unless ``aizynthfinder_stocks`` selects a subset.

    Returns:
        Path to the stock file.
    """
    from hedgehog.setup._aizynthfinder import ensure_aizynthfinder

    config_yml = ensure_aizynthfinder(project_root)
    stock_path = get_bld_stock_path(project_root)
    manifest = _manifest_path(stock_path)

    if refresh or not (stock_path.exists() and manifest.exists()):
        if not confirm_download(
            "BLDpharm stock", "110 MB from PubChem, 40 MB stock file"
        ):
            raise RuntimeError("BLDpharm stock download declined by user.")
        if refresh:
            shutil.rmtree(stock_path.parent, ignore_errors=True)
        stock_path.parent.mkdir(parents=True, exist_ok=True)
        build_bld_stock(stock_path, jobs=jobs)
    else:
        logger.info("BLDpharm stock found at %s", stock_path)

    register_stock(config_yml, BLD_STOCK_KEY, stock_path)
    return stock_path
