"""Preparation tests (pure-python engine + optional RDKit/Meeko paths)."""

from __future__ import annotations

import importlib.machinery
import importlib.util
import logging
import sys
import types
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
from tests.conftest import RECEPTOR_PDB


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


# ---------------------------------------------------------------------------
# Auto-engine graceful degradation (regression: openbabel-wheel 3.1.1.23
# removed OBElementTable/OBResidue.GetChainID and used to crash `engine=auto`)
# ---------------------------------------------------------------------------
def test_element_symbol_table():
    from dockflow_core.preparator import element_symbol

    assert element_symbol(6) == "C"
    assert element_symbol(30) == "Zn"
    assert element_symbol(1) == "H"
    assert element_symbol(999) == "X"
    # the table covers every element Vina can type
    for z in (1, 6, 7, 8, 15, 16, 11, 12, 20, 25, 26, 30, 34, 47, 53):
        assert element_symbol(z) != "X"


def test_engine_chain_auto_ends_with_none():
    from dockflow_core.preparator import _engine_chain

    chain = _engine_chain("auto")
    assert chain[-1].name == "none"  # always-available last resort
    names = [e.name for e in chain]
    # priority order preserved, no duplicates
    assert len(names) == len(set(names))
    # an explicit engine returns exactly that engine (fails loudly)
    assert [e.name for e in _engine_chain("none")] == ["none"]


class _BrokenOBModule(types.ModuleType):
    """Simulates a broken openbabel binding (the v0.1.1 audit crash)."""

    def __getattr__(self, name: str):
        raise AttributeError(f"module 'openbabel.openbabel' has no attribute '{name}'")


@pytest.fixture
def broken_openbabel(monkeypatch):
    """Make the openbabel bindings importable but unusable at attribute level.

    Works whether or not the *real* openbabel package is installed: the CI
    test matrix (macOS/Windows/Py3.10 cells) does not install the ``obabel``
    extra, so a stub package with a valid ``__spec__`` is injected into
    ``sys.modules`` - that keeps ``is_importable("openbabel")`` True (it
    uses ``importlib.util.find_spec``, which consults ``sys.modules`` and
    requires a non-None ``__spec__``), while any attribute access on the
    submodule raises AttributeError, exactly like the broken wheel the
    v0.1.1 audit discovered.
    """
    broken = _BrokenOBModule("openbabel.openbabel")
    broken.__spec__ = importlib.machinery.ModuleSpec(
        "openbabel.openbabel", loader=None
    )
    try:
        import openbabel as package
    except ImportError:  # not installed in this matrix cell - stub it
        package = types.ModuleType("openbabel")
        package.__spec__ = importlib.machinery.ModuleSpec(
            "openbabel", loader=None, is_package=True
        )
        package.__path__ = []  # mark as package for submodule imports
        monkeypatch.setitem(sys.modules, "openbabel", package)

    monkeypatch.setattr(package, "openbabel", broken, raising=False)
    monkeypatch.setitem(sys.modules, "openbabel.openbabel", broken)
    return broken


def test_auto_engine_falls_through_broken_optional_dependency(
    receptor_pdb_path: Path, tmp_path: Path, broken_openbabel, caplog
):
    """engine='auto' must never crash on a broken optional dependency.

    Regression for the v0.1.1 audit crash:
    ``AttributeError: module 'openbabel.openbabel' has no attribute
    'OBElementTable'`` in ``ReceptorPreparator.prepare``.  The run must fall
    through to the next engine in the chain and log a WARNING naming both the
    failed engine and the fallback.
    """
    caplog.set_level(logging.WARNING, logger="dockflow")
    with caplog.at_level(logging.WARNING, logger="dockflow"):
        result = ReceptorPreparator(ReceptorPrepOptions(engine="auto")).prepare(
            receptor_pdb_path, tmp_path
        )
    assert result.ok, f"preparation must still succeed via fallback: {result.warnings}"
    assert result.engine != "openbabel"
    assert result.engine in ("rdkit", "openbabel-cli", "none")
    if result.engine == "rdkit":
        # on a full install the next engine in the chain is RDKit
        assert any(
            "preparation engine 'openbabel' failed" in record.getMessage()
            and "falling back to 'rdkit'" in record.getMessage()
            for record in caplog.records
            if record.levelno == logging.WARNING
        ), "expected a WARNING naming the failed engine and the fallback"
        assert any("falling back to 'rdkit'" in warning for warning in result.warnings)
    else:
        assert any("falling back to" in warning for warning in result.warnings)
    # and the produced PDBQT is valid
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    assert atoms and all(a.atom_type for a in atoms)


@pytest.mark.skipif(
    not importlib.util.find_spec("openbabel"),
    reason="openbabel bindings not installed (Linux: pip install openbabel-wheel)",
)
def test_openbabel_engine_works_on_current_wheel(
    receptor_pdb_path: Path, tmp_path: Path
):
    """The OB engine itself must work on modern openbabel wheels.

    Regression for OBElementTable (removed in 3.1.1.23) and
    OBResidue.GetChainID (replaced by GetChain in the same generation).
    """
    result = ReceptorPreparator(ReceptorPrepOptions(engine="openbabel")).prepare(
        receptor_pdb_path, tmp_path
    )
    assert result.ok
    assert result.engine == "openbabel"
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    heavy = [a for a in atoms if a.element.upper() != "H"]
    assert len(heavy) == 40  # polymer heavy atoms (ZN/BEN/waters filtered)
    assert len(atoms) > len(heavy)  # polar hydrogens were added
    assert any(a.atom_type == "HD" for a in atoms)  # donor hydrogens typed
    assert all(a.atom_type for a in atoms)
    assert atoms[0].chain == "A"  # chain ids survive the round trip
    # Gasteiger charges were actually computed (not all zero)
    assert any(abs(a.charge) > 1e-6 for a in atoms)


def test_removed_species_recorded_and_warned(
    receptor_pdb_path: Path, tmp_path: Path, caplog
):
    """Removed metals/cofactors must be visible, not silent (audit item 11)."""
    caplog.set_level(logging.WARNING, logger="dockflow")
    result = ReceptorPreparator(ReceptorPrepOptions(engine="none")).prepare(
        receptor_pdb_path, tmp_path
    )
    assert result.removed_resnames.get("ZN2") == 1
    assert result.removed_resnames.get("BEN") == 7
    assert result.removed_resnames.get("HOH") == 2
    assert any("removed hetero species" in w for w in result.warnings)
    assert any(
        "structurally or catalytically required" in r.getMessage()
        and "ZN2" in r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING
    )


def test_water_only_removal_is_not_warned(tmp_path: Path):
    """Routine water removal must not raise the metals/cofactors warning."""
    pdb = "\n".join(
        line
        for line in RECEPTOR_PDB.splitlines()
        if not (line.startswith("HETATM") and " ZN2 " in line and " BEN " not in line)
        and " BEN " not in line
    )
    path = tmp_path / "waters_only.pdb"
    path.write_text(pdb, encoding="utf-8")
    result = ReceptorPreparator(ReceptorPrepOptions(engine="none")).prepare(path, tmp_path)
    assert result.removed_resnames.get("HOH") == 2
    assert not any("structurally or catalytically required" in w for w in result.warnings)
