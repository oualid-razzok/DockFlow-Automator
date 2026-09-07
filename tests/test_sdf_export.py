"""SDF pose export tests (audit item 29)."""

from __future__ import annotations

from pathlib import Path

import pytest

from dockflow_core.models import DockingResult, LigandRecord, PoseRecord
from dockflow_core.sdf_export import export_poses_sdf

try:
    from rdkit import Chem  # noqa: F401

    _HAS_RDKIT = True
except ImportError:  # pragma: no cover
    _HAS_RDKIT = False


def _result(tmp_path: Path, docked_pdbqt_text: str) -> DockingResult:
    out = tmp_path / "lig_out.pdbqt"
    out.write_text(docked_pdbqt_text, encoding="utf-8")
    return DockingResult(
        ligand=LigandRecord(identifier="lig", num_heavy_atoms=2),
        ligand_name="lig",
        poses=[
            PoseRecord(model=1, affinity=-9.423, rmsd_lb=0.0, rmsd_ub=0.0,
                       crystal_rmsd=1.25),
            PoseRecord(model=2, affinity=-8.711, rmsd_lb=1.2, rmsd_ub=2.1),
        ],
        out_path=out,
        backend="python",
    )


@pytest.mark.skipif(not _HAS_RDKIT, reason="rdkit not installed")
def test_export_poses_sdf_writes_one_record_per_pose(tmp_path: Path,
                                                     docked_pdbqt_text: str):
    result = _result(tmp_path, docked_pdbqt_text)
    target = export_poses_sdf(result, tmp_path)
    assert target is not None and target.name == "lig_poses.sdf"
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(target), removeHs=False)
    mols = [m for m in supplier if m is not None]
    assert len(mols) == 3  # the fixture has 3 MODELs
    # SD fields: scores + crystal_rmsd where known
    first = mols[0]
    assert first.GetProp("vina_affinity") == "-9.423"
    assert first.GetProp("vina_rmsd_lb") == "0.000"
    assert first.GetProp("crystal_rmsd") == "1.250"
    assert first.GetProp("docking_backend") == "python"
    assert first.GetProp("_Name") == "lig_pose1"
    # pose 2 has no crystal_rmsd -> field absent
    assert not mols[1].HasProp("crystal_rmsd")


@pytest.mark.skipif(not _HAS_RDKIT, reason="rdkit not installed")
def test_export_poses_sdf_missing_file(tmp_path: Path):
    result = DockingResult(ligand_name="ghost", out_path=tmp_path / "nope.pdbqt")
    assert export_poses_sdf(result, tmp_path) is None


def test_pose_pdb_block_conversion(tmp_path: Path, ligand_pdbqt_text: str):
    from dockflow_core.sdf_export import _pose_pdb_block

    source = tmp_path / "lig.pdbqt"
    source.write_text(ligand_pdbqt_text, encoding="utf-8")
    block = _pose_pdb_block(source)
    assert block.startswith("ATOM")
    assert "END" in block
    assert "ROOT" not in block and "TORSDOF" not in block
    assert "REMARK" not in block
