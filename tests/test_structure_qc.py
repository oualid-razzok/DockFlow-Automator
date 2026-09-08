"""Tests for receptor structure quality checks (work order items 6, 7, 10)
and engine comparability (peer item 4)."""

from dockflow_core.pdbio import Atom
from dockflow_core.preparator import engine_comparability
from dockflow_core.structure_qc import (
    detect_disulfides,
    detect_metal_coordination,
    detect_missing_residues,
    detect_reduced_cysteines,
    disulfide_partner_map,
    metal_warning,
    missing_residue_warning,
    reduced_cys_warning,
)


def _atom(name, resname, chain, resseq, x, y, z, element, record="ATOM"):
    return Atom(name=name, resname=resname, chain=chain, resseq=resseq,
                x=x, y=y, z=z, element=element, record_type=record,
                occupancy=1.0, bfactor=0.0)


def _backbone(resname, chain, resseq, x=0.0):
    """A minimal polymer residue (CA atom only) for gap walking."""
    return _atom("CA", resname, chain, resseq, x, 0.0, 0.0, "C")


# -- item 6: missing residues ------------------------------------------------
def test_missing_residues_detected():
    atoms = ([_backbone("ALA", "A", i) for i in range(1, 45)]
             + [_backbone("GLY", "A", i) for i in range(48, 52)])
    gaps = detect_missing_residues(atoms)
    assert gaps == [{"chain": "A", "resseq_range": "45-47", "count": 3}]


def test_missing_residues_insertion_codes_are_contiguous():
    # residue 11 with insertion code A directly follows 11 with '': no gap
    with_icode = [_backbone("ALA", "A", 10)]
    with_icode.append(Atom(name="CA", resname="GLY", chain="A", resseq=11,
                           icode="A", x=0, y=0, z=0, element="C",
                           record_type="ATOM", occupancy=1.0, bfactor=0.0))
    assert detect_missing_residues(with_icode) == []
    # numbering 10 -> 12 with nothing in between: one missing residue
    atoms = [_backbone("ALA", "A", 10), _backbone("GLY", "A", 12)]
    assert detect_missing_residues(atoms) == [
        {"chain": "A", "resseq_range": "11", "count": 1}]


def test_missing_residues_hetatoms_do_not_count():
    atoms = [_backbone("ALA", "A", 1), _backbone("GLY", "A", 2)]
    atoms.append(_atom("XK2", "XK2", "A", 5, 0, 0, 0, "C", record="HETATM"))
    assert detect_missing_residues(atoms) == []


def test_missing_residue_warning_text():
    warning = missing_residue_warning(
        [{"chain": "A", "resseq_range": "45-47", "count": 3}])
    assert "missing residues" in warning
    assert "chain A: 45-47" in warning
    assert "Modeller" in warning or "SWISS-MODEL" in warning
    assert missing_residue_warning([]) is None


# -- item 7: disulfides -------------------------------------------------------
def test_disulfide_detected():
    atoms = [
        _atom("SG", "CYS", "A", 14, 0.0, 0.0, 0.0, "S"),
        _atom("SG", "CYS", "A", 32, 2.0, 0.0, 0.0, "S"),
    ]
    bridges = detect_disulfides(atoms)
    assert bridges == [{"chain_a": "A", "resseq_a": 14,
                        "chain_b": "A", "resseq_b": 32,
                        "distance_a": 2.0}]


def test_disulfide_not_detected_when_far():
    atoms = [
        _atom("SG", "CYS", "A", 14, 0.0, 0.0, 0.0, "S"),
        _atom("SG", "CYS", "A", 32, 6.0, 0.0, 0.0, "S"),
    ]
    assert detect_disulfides(atoms) == []


def test_reduced_cysteines_flagged_near_partner():
    atoms = [
        _atom("SG", "CYS", "A", 145, 0.0, 0.0, 0.0, "S"),
        _atom("SG", "CYS", "B", 87, 5.0, 0.0, 0.0, "S"),
    ]
    flagged = detect_reduced_cysteines(atoms, [])
    assert flagged and flagged[0]["resseq"] == 145
    warning = reduced_cys_warning(flagged)
    assert "Cys145" in warning and "reduced" in warning


def test_disulfide_partner_map_symmetric():
    partners = disulfide_partner_map(
        [{"chain_a": "A", "resseq_a": 14, "chain_b": "B", "resseq_b": 87}])
    assert partners[("A", 14)] == ("B", 87)
    assert partners[("B", 87)] == ("A", 14)


# -- item 10: metals ----------------------------------------------------------
def test_metal_coordination_octahedral():
    atoms = [_atom("ZN", "ZN", "A", 401, 0.0, 0.0, 0.0, "ZN",
                   record="HETATM")]
    ligands = [
        ("HIS", 96, 2.0, 0.0, 0.0), ("HIS", 98, -2.0, 0.0, 0.0),
        ("ASP", 120, 0.0, 2.0, 0.0), ("ASP", 122, 0.0, -2.0, 0.0),
        ("H2O", 501, 0.0, 0.0, 2.0), ("H2O", 502, 0.0, 0.0, -2.0),
    ]
    for i, (resname, resseq, x, y, z) in enumerate(ligands):
        element = "O" if resname == "H2O" else "N"
        record = "HETATM" if resname == "H2O" else "ATOM"
        atoms.append(_atom(f"X{i}", resname, "A", resseq, x, y, z, element,
                           record=record))
    findings = detect_metal_coordination(atoms)
    assert len(findings) == 1
    finding = findings[0]
    assert finding["metal_resname"] == "ZN"
    assert finding["coordination_number"] == 4  # waters excluded
    assert finding["geometry_heuristic"] == "tetrahedral"
    residues = {(r["resname"], r["resseq"]) for r in finding["coordinating_residues"]}
    assert ("HIS", 96) in residues and ("ASP", 120) in residues


def test_metal_warning_mentions_scoring_limitation():
    findings = [{"metal_resname": "ZN", "geometry_heuristic": "tetrahedral",
                 "coordination_number": 4, "coordinating_residues": []}]
    warning = metal_warning(findings, backend="vina")
    assert "Vina has no metal-coordination scoring term" in warning
    assert "GNINA" in warning
    gnina_warning = metal_warning(findings, backend="gnina")
    assert "no metal-coordination scoring term" not in gnina_warning
    assert metal_warning([], backend="vina") is None


# -- item 4: engine comparability --------------------------------------------
def test_engine_comparability_rdkit():
    info = engine_comparability("rdkit")
    assert info["resolved_engine"] == "rdkit"
    assert info["comparable_to"] == ["rdkit"]
    assert "openbabel" in info["incomparable_to"]
    assert "none" in info["incomparable_to"]


def test_engine_comparability_none_incomparable_to_everything():
    info = engine_comparability("none (hydrogens disabled)")
    assert info["resolved_engine"] == "none"
    assert info["comparable_to"] == []
    assert info["incomparable_to"] == ["openbabel", "openbabel-cli", "rdkit"]


def test_engine_comparability_openbabel_family():
    info = engine_comparability("openbabel-cli")
    assert set(info["comparable_to"]) == {"openbabel", "openbabel-cli"}
    assert "rdkit" in info["incomparable_to"]
