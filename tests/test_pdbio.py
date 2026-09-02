"""PDB/PDBQT parsing and writing tests."""

from __future__ import annotations

from pathlib import Path

from dockflow_core.pdbio import (
    Atom,
    element_from_atom_name,
    filter_atoms,
    format_pdb_line,
    format_pdbqt_line,
    het_resnames,
    parse_pdb,
    parse_pdb_models,
    parse_pdbqt,
    parse_pdbqt_results,
    split_pdbqt_models,
    write_pdb,
    write_pdbqt,
)


def test_parse_receptor_basic(receptor_pdb_text: str):
    atoms = parse_pdb(receptor_pdb_text)
    # 40 polymer + ZN + 2 waters + benzene (6 C + 2 O1 altloc entries) = 51
    assert len(atoms) == 51
    first = atoms[0]
    assert first.name == "N" and first.resname == "ALA"
    assert first.chain == "A" and first.resseq == 1
    assert first.element == "N"
    zn = atoms[40]
    assert zn.resname == "ZN2" and zn.element == "Zn" and zn.record_type == "HETATM"
    ben = atoms[43]
    assert ben.resname == "BEN"


def test_parse_multi_model_pdb():
    text = (
        "MODEL     1\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ENDMDL\n"
        "MODEL     2\n"
        "ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N\n"
        "ENDMDL\n"
    )
    models = parse_pdb_models(text)
    assert len(models) == 2
    assert models[1][0].x == 1.0


def test_pdbqt_line_column_layout():
    atom = Atom(
        serial=1, name="N", resname="VAL", chain="A", resseq=1,
        x=13.199, y=24.047, z=45.489, occupancy=1.0, bfactor=0.0,
        element="N", charge=0.343, atom_type="N",
    )
    line = format_pdbqt_line(atom)
    assert len(line) == 79
    assert line[0:6] == "ATOM  "
    assert line[12:16] == " N  "
    assert line[17:20] == "VAL"
    assert line[21:22] == "A"
    assert line[22:26] == "   1"
    assert line[70:76].strip() == "0.343"
    assert line[77:79].strip() == "N"


def test_pdbqt_line_blank_chain_altloc():
    atom = Atom(serial=2, name="CA", resname="LIG", chain="", resseq=2,
                x=14.0, y=25.0, z=46.0, element="C", charge=-0.123, atom_type="A")
    line = format_pdbqt_line(atom)
    assert len(line) == 79
    assert line[16:17] == " "  # altloc column exists
    assert line[21:22] == " "  # chain column exists
    assert line[17:20] == "LIG"


def test_pdbqt_roundtrip(receptor_pdb_path: Path, tmp_path: Path):
    atoms = parse_pdb(receptor_pdb_path)
    typed = [
        Atom(
            serial=i + 1, name=a.name, resname=a.resname, chain=a.chain,
            resseq=a.resseq, x=a.x, y=a.y, z=a.z, element=a.element,
            occupancy=a.occupancy, bfactor=a.bfactor,
            charge=0.15 if i % 2 else -0.12,
            atom_type="C" if a.element == "C" else a.element,
            record_type=a.record_type,
        )
        for i, a in enumerate(atoms)
    ]
    out = tmp_path / "roundtrip.pdbqt"
    write_pdbqt(typed, out)
    parsed = parse_pdbqt(out).atoms
    assert len(parsed) == len(typed)
    for original, again in zip(typed, parsed, strict=True):
        assert again.name == original.name
        assert again.resname == original.resname
        assert again.chain == original.chain
        assert again.resseq == original.resseq
        assert abs(again.x - original.x) < 1e-6
        assert abs(again.charge - original.charge) < 1e-6
        assert again.atom_type == original.atom_type


def test_write_and_parse_pdb(receptor_pdb_path: Path, tmp_path: Path):
    atoms = parse_pdb(receptor_pdb_path)
    out = tmp_path / "written.pdb"
    write_pdb(atoms, out)
    again = parse_pdb(out)
    assert len(again) == len(atoms)
    assert again[0].resname == "ALA"
    assert again[0].chain == "A"


def test_parse_docked_pdbqt_poses(docked_pdbqt_text: str):
    data = parse_pdbqt(docked_pdbqt_text)
    assert len(data.models) == 3
    results = data.vina_results()
    assert [r.affinity for r in results] == [-9.423, -8.711, -7.905]
    assert results[1].rmsd_lb == 1.234 and results[1].rmsd_ub == 2.100
    poses = parse_pdbqt_results(docked_pdbqt_text)
    assert len(poses) == 3
    assert poses[0].affinity == -9.423


def test_split_pdbqt_models(docked_pdbqt_path: Path, tmp_path: Path):
    files = split_pdbqt_models(docked_pdbqt_path, tmp_path / "poses", "lig")
    assert len(files) == 3
    first = parse_pdbqt(files[0])
    assert first.models[0].vina_result is not None
    assert first.models[0].vina_result.affinity == -9.423


def test_ligand_pdbqt_parse(ligand_pdbqt_text: str):
    data = parse_pdbqt(ligand_pdbqt_text)
    assert data.torsdof == 2
    assert len(data.atoms) == 5
    oa = data.atoms[1]
    assert oa.atom_type == "OA"
    assert abs(oa.charge - (-0.641)) < 1e-9


def test_filter_atoms(receptor_pdb_path: Path):
    atoms = parse_pdb(receptor_pdb_path)
    without_water = filter_atoms(atoms)
    assert all(not a.is_water for a in without_water)
    assert len(without_water) < len(atoms)
    chain_b = filter_atoms(atoms, chains=["B"])
    assert chain_b == []
    with_zn = filter_atoms(atoms, keep_resnames=["ZN2"])
    assert any(a.resname == "ZN2" for a in with_zn)
    polymers_only = filter_atoms(atoms, keep_hetero=False)
    assert all(a.is_polymer for a in polymers_only)


def test_het_resnames(receptor_pdb_path: Path):
    atoms = parse_pdb(receptor_pdb_path)
    codes = het_resnames(atoms)
    assert "ZN2" in codes and "BEN" in codes
    assert "HOH" not in codes
    assert "ALA" not in codes


def test_element_guessing():
    assert element_from_atom_name("N", "ALA") == "N"
    assert element_from_atom_name("CA", "ALA") == "C"
    assert element_from_atom_name("CA", "ZN2") == "Ca"
    assert element_from_atom_name("ZN", "ZN2") == "Zn"
    assert element_from_atom_name("CL", "LIG") == "Cl"
    assert element_from_atom_name("HD11", "LEU") == "H"


def test_atom_name_field_alignment():
    # short names shift right by one column (" N  ")
    line = format_pdb_line(Atom(name="N", resname="VAL", chain="A", resseq=1))
    assert line[12:16] == " N  "
    # 4-character names fill the field completely
    line = format_pdb_line(Atom(name="HD11", resname="LEU", chain="A", resseq=2))
    assert line[12:16] == "HD11"
