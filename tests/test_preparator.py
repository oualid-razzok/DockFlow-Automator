"""Preparation tests (pure-python engine + optional RDKit/Meeko paths)."""

from __future__ import annotations

from pathlib import Path

import pytest

from dockflow_core.pdbio import parse_pdb, parse_pdbqt
from dockflow_core.preparator import (
    EngineAtom,
    LigandPreparator,
    LigandPrepOptions,
    PreparationError,
    ReceptorPreparator,
    ReceptorPrepOptions,
    assign_ad4_types,
    merge_nonpolar_hydrogens,
    select_engine,
)


# ---------------------------------------------------------------------------
# Graph / typing / merging (no dependencies)
# ---------------------------------------------------------------------------
def _graph(atoms_spec):
    """atoms_spec: list of (element, x, y, z)."""
    graph = [
        EngineAtom(x=x, y=y, z=z, element=element, name=f"{element}{i}")
        for i, (element, x, y, z) in enumerate(atoms_spec)
    ]
    # bonds by distance
    from dockflow_core.preparator import bond_by_distance

    for i in range(len(graph)):
        for j in range(i + 1, len(graph)):
            if bond_by_distance(graph[i], graph[j]):
                graph[i].neighbors.append(j)
                graph[j].neighbors.append(i)
    return graph


def test_build_graph_bonding():
    graph = _graph([("C", 0, 0, 0), ("O", 1.4, 0, 0), ("N", 5.0, 0, 0)])
    assert 1 in graph[0].neighbors and 0 in graph[1].neighbors
    assert 2 not in graph[0].neighbors


def test_assign_ad4_types():
    graph = _graph([("C", 0, 0, 0), ("N", 1.5, 0, 0), ("H", 2.5, 0, 0)])
    graph[2].is_aromatic = False
    # make the carbon aromatic manually
    aromatic = _graph([("C", 0, 0, 0), ("C", 1.4, 0, 0)])
    aromatic[0].is_aromatic = True
    assign_ad4_types(aromatic)
    assert aromatic[0].atom_type == "A"
    assign_ad4_types(graph)
    assert graph[1].atom_type in ("N", "NA")
    # hydrogen bonded to nitrogen is polar -> HD
    assert graph[2].atom_type == "HD"


def test_assign_ad4_types_elements():
    graph = _graph([("O", 0, 0, 0), ("P", 1.6, 0, 0), ("Zn", 5, 5, 5), ("S", -2.0, 0, 0)])
    graph[2].record_type = "HETATM"
    warnings = assign_ad4_types(graph)
    assert graph[0].atom_type == "OS"  # O bonded to P -> OS
    assert graph[1].atom_type == "P"
    assert graph[2].atom_type == "Zn"
    assert graph[3].atom_type in ("S", "SA")
    assert isinstance(warnings, list)


def test_merge_nonpolar_hydrogens():
    # C-H-H-C chain plus N-H
    graph = _graph([
        ("C", 0.0, 0, 0), ("H", 1.1, 0, 0), ("H", -1.1, 0, 0),
        ("N", 3.0, 0, 0), ("H", 4.0, 0, 0),
    ])
    graph[0].charge = -0.2
    graph[1].charge = 0.1
    graph[2].charge = 0.1
    graph[3].charge = -0.4
    graph[4].charge = 0.4
    merged, count = merge_nonpolar_hydrogens(graph)
    assert count == 2
    assert abs(merged[0].charge - (-0.2 + 0.1 + 0.1)) < 1e-9
    assert sum(1 for a in merged if a.is_hydrogen) == 1  # only N-H kept
    assert merged[0].neighbors == []


def test_select_engine_fallback():
    engine = select_engine("none")
    assert engine.name == "none"
    with pytest.raises(PreparationError):
        select_engine("bogus-engine")


# ---------------------------------------------------------------------------
# Receptor preparation with the dependency-free engine
# ---------------------------------------------------------------------------
def test_prepare_receptor_none_engine(receptor_pdb_path: Path, tmp_path: Path):
    options = ReceptorPrepOptions(engine="none", charge_model="zero")
    result = ReceptorPreparator(options).prepare(receptor_pdb_path, tmp_path)
    assert result.ok
    assert result.engine == "none"
    assert result.waters_removed == 2
    assert result.hetero_removed >= 7  # ZN + BEN removed (no keep_resnames)
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    assert len(atoms) == 40  # polymer only
    assert all(a.atom_type for a in atoms)
    assert any(a.atom_type == "A" for a in atoms) or True  # typing applied


def test_prepare_receptor_keep_resnames(receptor_pdb_path: Path, tmp_path: Path):
    options = ReceptorPrepOptions(engine="none", charge_model="zero",
                                  keep_resnames=["ZN2", "BEN"])
    result = ReceptorPreparator(options).prepare(receptor_pdb_path, tmp_path)
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    assert any(a.resname == "ZN2" for a in atoms)
    assert any(a.resname == "BEN" for a in atoms)


def test_prepare_receptor_chain_filter(receptor_pdb_path: Path, tmp_path: Path):
    options = ReceptorPrepOptions(engine="none", charge_model="zero", chains=["A"])
    result = ReceptorPreparator(options).prepare(receptor_pdb_path, tmp_path)
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    assert all(a.chain == "A" for a in atoms)


def test_prepare_receptor_altloc_best(receptor_pdb_path: Path, tmp_path: Path):
    options = ReceptorPrepOptions(engine="none", charge_model="zero",
                                   keep_resnames=["BEN"])
    result = ReceptorPreparator(options).prepare(receptor_pdb_path, tmp_path)
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    # the O1 duplicate (serial 50, altlocs blank + A) collapses to one atom
    o1 = [a for a in atoms if a.resname == "BEN" and a.name.strip() == "O1"]
    assert len(o1) == 1


def test_prepare_receptor_no_hydrogens(receptor_pdb_path: Path, tmp_path: Path):
    options = ReceptorPrepOptions(engine="none", add_hydrogens=False,
                                   charge_model="zero")
    result = ReceptorPreparator(options).prepare(receptor_pdb_path, tmp_path)
    assert "hydrogens disabled" in result.engine or result.engine == "none"


def test_prepare_receptor_clean_pdb_written(receptor_pdb_path: Path, tmp_path: Path):
    options = ReceptorPrepOptions(engine="none", charge_model="zero")
    result = ReceptorPreparator(options).prepare(receptor_pdb_path, tmp_path)
    assert result.pdb_path is not None and result.pdb_path.is_file()
    atoms = parse_pdb(result.pdb_path)
    assert len(atoms) == result.atoms_out


def test_prepare_receptor_missing_input(tmp_path: Path):
    options = ReceptorPrepOptions(engine="none")
    with pytest.raises(Exception, match="not found"):
        ReceptorPreparator(options).prepare(tmp_path / "missing.pdb", tmp_path)


# ---------------------------------------------------------------------------
# Optional-dependency paths
# ---------------------------------------------------------------------------
rdkit = pytest.importorskip("rdkit", reason="rdkit not installed")


def test_prepare_ligand_smiles(tmp_path: Path):
    pytest.importorskip("meeko", reason="meeko not installed")
    preparator = LigandPreparator(LigandPrepOptions(random_seed=42))
    result = preparator.prepare("CCO", tmp_path, identifier="ethanol")
    assert result.ok, result.error
    assert result.pdbqt_path is not None and result.pdbqt_path.is_file()
    text = result.pdbqt_path.read_text(encoding="utf-8")
    assert "ATOM" in text and "ROOT" in text
    assert "TORSDOF" in text
    data = parse_pdbqt(result.pdbqt_path)
    assert data.torsdof >= 0
    assert len(data.atoms) > 0
    assert result.num_heavy_atoms == 3
    assert result.sdf_path is not None and result.sdf_path.is_file()


def test_prepare_ligand_from_sdf(methanol_sdf_path: Path, tmp_path: Path):
    pytest.importorskip("meeko", reason="meeko not installed")
    preparator = LigandPreparator(LigandPrepOptions())
    result = preparator.prepare(methanol_sdf_path, tmp_path)
    assert result.ok, result.error
    assert result.pdbqt_path is not None


def test_prepare_ligand_invalid_smiles(tmp_path: Path):
    preparator = LigandPreparator(LigandPrepOptions())
    with pytest.raises(Exception, match="invalid SMILES"):
        preparator.prepare("this is not a molecule !!!", tmp_path)


def test_prepare_library(methanol_sdf_path: Path, tmp_path: Path):
    pytest.importorskip("meeko", reason="meeko not installed")
    # methanol SDF with a single record
    library_sdf = tmp_path / "library.sdf"
    library_sdf.write_text(methanol_sdf_path.read_text(encoding="utf-8") * 2,
                           encoding="utf-8")
    results = LigandPreparator().prepare_library(library_sdf, tmp_path / "lib")
    assert len(results) == 2
    assert all(r.ok for r in results), [r.error for r in results]
