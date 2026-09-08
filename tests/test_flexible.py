"""Tests for flexible receptor residues (peer item 8)."""

from pathlib import Path

import pytest

from dockflow_core.flexible import (
    CHI_BONDS,
    auto_flexible_residues,
    build_flexible_pdbqts,
    parse_flexible_residues,
)
from dockflow_core.utils import DockFlowError


def _pdbqt_atom(serial, name, resname, chain, resseq, x, y, z, element="C",
                atom_type=None):
    atom_type = atom_type or element
    return (f"ATOM  {serial:5d} {name:<4s} {resname:>3s} {chain:1s}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    0.100 {atom_type:<2s}")


def _receptor(tmp_path: Path) -> Path:
    """ASP25 (side chain), LYS30 (long chain), GLY40 (no side chain)."""
    lines = []
    serial = 1
    # ASP 25: N, CA, C, O backbone + CB, CG, OD1, OD2
    positions = {
        "N": (0.0, 0.0, 0.0), "CA": (1.5, 0.0, 0.0), "C": (2.2, 1.3, 0.0),
        "O": (1.6, 2.4, 0.0), "CB": (2.3, -1.2, 0.0), "CG": (3.8, -1.4, 0.0),
        "OD1": (4.5, -0.4, 0.0), "OD2": (4.3, -2.6, 0.0),
    }
    for name, (x, y, z) in positions.items():
        element = "O" if name.startswith(("O", "OD")) else ("N" if name == "N" else "C")
        at = "OA" if name.startswith("OD") else element
        lines.append(_pdbqt_atom(serial, name, "ASP", "A", 25, x, y, z,
                                 element, at))
        serial += 1
    # LYS 30 with a flexible chain
    lys = {"N": (10, 0, 0), "CA": (11.5, 0, 0), "C": (12.2, 1.3, 0),
           "O": (11.6, 2.4, 0), "CB": (12.3, -1.2, 0), "CG": (13.8, -1.4, 0),
           "CD": (14.5, -2.7, 0), "CE": (16.0, -2.9, 0), "NZ": (16.7, -4.2, 0)}
    for name, (x, y, z) in lys.items():
        element = "N" if name in ("N", "NZ") else ("O" if name == "O" else "C")
        at = "NA" if name == "NZ" else element
        lines.append(_pdbqt_atom(serial, name, "LYS", "A", 30, x, y, z,
                                 element, at))
        serial += 1
    # GLY 40: no side chain
    for name, (x, y, z) in {"N": (30, 0, 0), "CA": (31.5, 0, 0),
                            "C": (32.2, 1.3, 0), "O": (31.6, 2.4, 0)}.items():
        element = "N" if name == "N" else ("O" if name == "O" else "C")
        lines.append(_pdbqt_atom(serial, name, "GLY", "A", 40, x, y, z, element))
        serial += 1
    path = tmp_path / "receptor.pdbqt"
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


def test_build_flexible_pdbqts_asp_lys(tmp_path):
    rigid, flex, report = build_flexible_pdbqts(
        _receptor(tmp_path),
        [{"chain": "A", "resseq": 25}, {"chain": "A", "resseq": 30}],
        tmp_path)
    assert rigid.is_file() and flex.is_file()
    flex_text = flex.read_text()
    assert "BEGIN_RES ASP A 25" in flex_text
    assert "BEGIN_RES LYS A 30" in flex_text
    # root is CB; BRANCH lines reference serials
    assert flex_text.splitlines()[1] == "ROOT"
    branch_lines = [line for line in flex_text.splitlines()
                    if line.startswith("BRANCH")]
    end_branch = [line for line in flex_text.splitlines()
                  if line.startswith("ENDBRANCH")]
    assert len(branch_lines) == 1 + len(CHI_BONDS["LYS"])  # ASP chi2 + LYS chi1..3
    assert len(branch_lines) == len(end_branch)
    for line in branch_lines:
        parts = line.split()
        assert parts[1].isdigit() and parts[2].isdigit()
    # rigid part keeps backbone, drops the flexible side chains
    rigid_text = rigid.read_text()
    assert " CA  ASP A  25" in rigid_text
    assert " CB  ASP A  25" not in rigid_text
    assert " CG  LYS A  30" not in rigid_text
    assert " CA  LYS A  30" in rigid_text
    flexible = [r for r in report if r["status"] == "flexible"]
    assert {r["resname"] for r in flexible} == {"ASP", "LYS"}


def test_build_flexible_skips_glycine_with_reason(tmp_path):
    rigid, flex, report = build_flexible_pdbqts(
        _receptor(tmp_path), [{"chain": "A", "resseq": 40},
                              {"chain": "A", "resseq": 25}], tmp_path)
    skipped = [r for r in report if r["resseq"] == 40]
    assert skipped and skipped[0]["status"] == "skipped"
    assert "no rotatable" in skipped[0]["reason"]
    assert flex.is_file()  # ASP still flexible


def test_build_flexible_residue_not_found(tmp_path):
    with pytest.raises(DockFlowError, match="no flexible side chains"):
        build_flexible_pdbqts(_receptor(tmp_path),
                              [{"chain": "Z", "resseq": 999}], tmp_path)


def test_parse_flexible_residues_forms():
    assert parse_flexible_residues(None) == []
    assert parse_flexible_residues([]) == []
    assert parse_flexible_residues([{"chain": "A", "resseq": 25}]) == [
        {"chain": "A", "resseq": 25}]
    auto = parse_flexible_residues("auto")
    assert auto == [{"auto": True}]
    with pytest.raises(DockFlowError):
        parse_flexible_residues("banana")
    with pytest.raises(DockFlowError):
        parse_flexible_residues([{"chain": "A"}])  # resseq missing


def test_auto_flexible_residues_contacts(tmp_path):
    receptor = _receptor(tmp_path)
    ligand = tmp_path / "lig.pdb"
    ligand.write_text(
        "ATOM      1  C   XK2 A 500      4.300  -1.400   0.000"
        "  1.00  0.00           C\n", encoding="utf-8")
    residues = auto_flexible_residues(receptor, ligand, contact_cutoff=4.5)
    chains = {(r["chain"], r["resseq"]) for r in residues}
    assert ("A", 25) in chains   # ligand sits on the ASP side chain
    assert ("A", 40) not in chains  # GLY far away


def test_flex_report_records_rotatable_bond_counts(tmp_path):
    _, _, report = build_flexible_pdbqts(
        _receptor(tmp_path), [{"chain": "A", "resseq": 30}], tmp_path)
    lys = [r for r in report if r["resseq"] == 30][0]
    assert lys["rotatable_bonds"] == len(CHI_BONDS["LYS"])
    assert lys["atoms_flexible"] == 5  # CB, CG, CD, CE, NZ
