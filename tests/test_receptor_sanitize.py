"""Tests for receptor PDB sanitization for GNINA."""

from pathlib import Path

from hedgehog.docking.receptor_sanitize import (
    ensure_gnina_compatible_receptor,
    sanitize_receptor_pdb,
)


MAESTRO_SAMPLE = """\
ANISOU  839  NZ ALYS A  56     3449   1983   3690   -789  -1310   -278
ATOM    839  NZ ALYS A  56      12.885   6.281  29.997  0.50 24.01           N1+
ATOM    840  NZ BLYS A  56      10.494   7.366  30.904  0.50 29.74           N1+
ATOM      1  N   MET A   1      14.802  -3.618  40.125  1.00 32.25           N1+
ATOM     28  OE2 GLU A   2      14.527  -1.779  42.138  0.50 34.76           O1-
"""


def test_sanitize_maestro_elements_and_alt_conformers(tmp_path: Path):
    source = tmp_path / "maestro.pdb"
    source.write_text(MAESTRO_SAMPLE, encoding="utf-8")
    destination = tmp_path / "clean.pdb"

    sanitize_receptor_pdb(source, destination)

    lines = destination.read_text(encoding="utf-8").splitlines()
    assert all(not line.startswith("ANISOU") for line in lines)
    assert len(lines) == 3

    nz_lines = [line for line in lines if "NZ" in line and "LYS" in line]
    assert len(nz_lines) == 1
    assert nz_lines[0][76:78] == " N"
    assert nz_lines[0][78:80] == "1+"

    met_n = next(line for line in lines if "MET" in line)
    assert met_n[76:78] == " N"
    assert met_n[78:80] == "1+"

    glu_o = next(line for line in lines if "OE2" in line)
    assert glu_o[76:78] == " O"
    assert glu_o[78:80] == "1-"


def test_ensure_gnina_compatible_receptor_caches_in_workdir(tmp_path: Path):
    source = tmp_path / "receptor.pdb"
    source.write_text(MAESTRO_SAMPLE, encoding="utf-8")
    workdir = tmp_path / "docking"

    first = ensure_gnina_compatible_receptor(source, workdir)
    second = ensure_gnina_compatible_receptor(source, workdir)

    assert first == second
    assert first.endswith("receptor_gnina.pdb")
    assert Path(first).exists()
