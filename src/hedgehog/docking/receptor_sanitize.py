"""Sanitize receptor PDB files for GNINA / AutoDock parsers.

Maestro-exported structures often use non-standard element notation (``N1+``,
``O1-``) in columns 77-80, interleaved ``ANISOU`` records, and alternate
conformations. GNINA's receptor parser treats these as invalid AutoDock atom
types (e.g. ``NA``); SMINA is more tolerant.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from hedgehog.configs.logger import logger

_MAESTRO_ELEMENT_MAP = {
    "N1+": ("N", "1+"),
    "O1-": ("O", "1-"),
}


@dataclass(frozen=True)
class _AtomRecord:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: str
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    charge: str

    @property
    def dedup_key(self) -> tuple[str, str, str, str]:
        return (self.chain, self.resseq, self.resname, self.name.strip())

    @property
    def altloc_rank(self) -> int:
        if not self.altloc:
            return 0
        if self.altloc == "A":
            return 1
        return 2


def _parse_atom_line(line: str) -> _AtomRecord | None:
    if len(line) < 66 or not line.startswith(("ATOM", "HETATM")):
        return None

    element_tail = line[76:].strip() if len(line) > 76 else ""
    element, charge = _normalize_element_and_charge(element_tail, line[12:16])

    return _AtomRecord(
        record=line[0:6].strip(),
        serial=int(line[6:11]),
        name=line[12:16],
        altloc=line[16] if line[16] != " " else "",
        resname=line[17:20].strip(),
        chain=line[21],
        resseq=line[22:26],
        icode=line[26] if len(line) > 26 and line[26] != " " else "",
        x=float(line[30:38]),
        y=float(line[38:46]),
        z=float(line[46:54]),
        occupancy=float(line[54:60]),
        bfactor=float(line[60:66]),
        element=element,
        charge=charge,
    )


def _normalize_element_and_charge(element_tail: str, atom_name: str) -> tuple[str, str]:
    if element_tail in _MAESTRO_ELEMENT_MAP:
        return _MAESTRO_ELEMENT_MAP[element_tail]

    if len(element_tail) == 2 and element_tail[1] in "+-":
        return element_tail[0], element_tail[1:]

    if element_tail:
        return element_tail[:2].rjust(2), ""

    guess = atom_name.strip()[0:1]
    return guess.rjust(2), ""


def _format_atom_line(atom: _AtomRecord, serial: int) -> str:
    record = atom.record.ljust(6)[:6]
    name = atom.name if len(atom.name) == 4 else atom.name.rjust(4)
    altloc = atom.altloc or " "
    resname = atom.resname.rjust(3)[:3]
    chain = atom.chain or " "
    resseq = atom.resseq.rjust(4)[:4]
    icode = atom.icode or " "
    element = atom.element.rjust(2)[:2]
    charge = atom.charge.ljust(2)[:2]

    return (
        f"{record}{serial:5d} {name}{altloc}{resname} {chain}{resseq}{icode}   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          "
        f"{element}{charge}"
    )


def _select_alt_conformers(atoms: list[_AtomRecord]) -> list[_AtomRecord]:
    best: dict[tuple[str, str, str, str], _AtomRecord] = {}
    for atom in atoms:
        key = atom.dedup_key
        current = best.get(key)
        if current is None:
            best[key] = atom
            continue
        if atom.occupancy > current.occupancy:
            best[key] = atom
            continue
        if atom.occupancy == current.occupancy and atom.altloc_rank < current.altloc_rank:
            best[key] = atom
    return list(best.values())


def sanitize_receptor_pdb(source: Path, destination: Path) -> Path:
    """Write a GNINA-compatible receptor PDB from a Maestro-style source file."""
    atoms: list[_AtomRecord] = []
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parsed = _parse_atom_line(line)
        if parsed is not None:
            atoms.append(parsed)

    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in receptor PDB: {source}")

    atoms = _select_alt_conformers(atoms)
    atoms.sort(key=lambda atom: (atom.chain, int(atom.resseq), atom.name, atom.altloc))

    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = [_format_atom_line(atom, serial) for serial, atom in enumerate(atoms, start=1)]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def ensure_gnina_compatible_receptor(
    receptor_path: str | Path,
    workdir: str | Path,
) -> str:
    """Return a GNINA-safe receptor path, sanitizing Maestro exports when needed."""
    source = Path(receptor_path)
    if not source.exists():
        return str(source)

    destination = Path(workdir) / "_workdir" / "receptor_gnina.pdb"
    if destination.exists() and destination.stat().st_mtime >= source.stat().st_mtime:
        return str(destination.resolve())

    try:
        sanitize_receptor_pdb(source, destination)
    except Exception as exc:
        logger.warning(
            "GNINA receptor sanitization failed for %s (%s); using original file",
            source,
            exc,
        )
        return str(source.resolve())

    logger.info(
        "GNINA: Sanitized receptor for AutoDock compatibility: %s -> %s",
        source,
        destination,
    )
    return str(destination.resolve())
