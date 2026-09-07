"""Export docked poses from PDBQT to SDF (audit item 29).

Vina writes multi-pose ``*_out.pdbqt`` files which most visualisation and
chemistry tooling cannot read directly.  This module converts them to
``<ligand>_poses.sdf`` with **one record per pose** and SD data fields for
the scores:

* ``vina_affinity``   - predicted binding affinity (kcal/mol)
* ``vina_rmsd_lb`` / ``vina_rmsd_ub`` - RMSD bounds to the best mode
* ``crystal_rmsd``    - symmetry-tolerant RMSD to the co-crystal pose
  (present only for redocking runs, filled by the analysis stage)

Two conversion strategies, most rigorous first:

1. **Meeko** (``PDBQTMolecule`` + ``RDKitMolCreate``): reconstructs the
   proper bond orders and stereochemistry recorded in the PDBQT remarks.
2. **RDKit PDB fallback**: per-pose ``MolFromPDBBlock`` with proximity
   bonding - used when Meeko is unavailable or rejects the file.

Enabled by default via ``docking.export_sdf: true``.
"""

from __future__ import annotations

from pathlib import Path

from .models import DockingResult
from .pdbio import split_pdbqt_models
from .utils import DockFlowError, get_logger

logger = get_logger("sdf_export")

__all__ = ["export_poses_sdf"]


def export_poses_sdf(result: DockingResult, out_dir: str | Path) -> Path | None:
    """Write ``<ligand>_poses.sdf`` with one record per docked pose.

    Returns the path, or ``None`` when the file could not be produced (a
    warning is logged; SDF export is a convenience, never a hard failure).
    """
    if result.out_path is None or not Path(result.out_path).is_file():
        return None
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{result.ligand_name}_poses.sdf"
    mols = _mols_via_meeko(result.out_path) or _mols_via_rdkit_fallback(result)
    if not mols:
        logger.warning("SDF export produced no molecules for %s", result.ligand_name)
        return None
    _write_sdf(mols, result, target)
    logger.info("exported %d pose(s) to %s", len(mols), target)
    return target


def _mols_via_meeko(pdbqt_path: Path):
    """Reconstruct molecules with correct bond orders via Meeko."""
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
    except ImportError:
        return []
    try:
        pdbqt_mol = PDBQTMolecule.from_file(str(pdbqt_path), skip_typing=True)
        return RDKitMolCreate.from_pdbqt_mol(pdbqt_mol) or []
    except Exception as exc:  # noqa: BLE001 - malformed input, fall back
        logger.debug("meeko SDF conversion failed (%s); using PDB fallback", exc)
        return []


def _mols_via_rdkit_fallback(result: DockingResult):
    """Per-pose RDKit molecules with proximity-bonded topology."""
    try:
        from rdkit import Chem
    except ImportError:
        return []
    if result.out_path is None:
        return []
    pose_files = split_pdbqt_models(result.out_path, Path(result.out_path).parent,
                                    result.ligand_name or "pose")
    mols = []
    for pose_file in pose_files:
        block = _pose_pdb_block(pose_file)
        if not block:
            continue
        mol = Chem.MolFromPDBBlock(block, removeHs=False, proximityBonding=True)
        if mol is not None:
            mols.append(mol)
    return mols


_PDBQT_ATOM_PREFIXES = ("ATOM", "HETATM")


def _pose_pdb_block(pose_file: Path) -> str:
    """Convert one single-pose PDBQT file back into a PDB block."""
    lines = []
    for line in pose_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(_PDBQT_ATOM_PREFIXES):
            # PDBQT columns 1-70 are plain PDB; drop charge/type columns.
            lines.append(line[:66] + line[76:78].rjust(2))
        elif line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF",
                              "REMARK", "MODEL", "ENDMDL", "BEGIN_RES", "END_RES")):
            continue
        else:
            lines.append(line)
    lines.append("END")
    return "\n".join(lines) + "\n"


def _write_sdf(mols, result: DockingResult, target: Path) -> None:
    """Write each molecule/conformer as one SDF record with SD score fields."""
    from rdkit import Chem

    writer = Chem.SDWriter(str(target))
    try:
        writer.SetKekulize(True)
        for index, mol in enumerate(mols):
            record = Chem.Mol(mol)
            poses = result.poses[index] if index < len(result.poses) else None
            record.SetProp("_Name", f"{result.ligand_name}_pose{index + 1}")
            if poses is not None:
                record.SetProp("vina_affinity", f"{poses.affinity:.3f}")
                record.SetProp("vina_rmsd_lb", f"{poses.rmsd_lb:.3f}")
                record.SetProp("vina_rmsd_ub", f"{poses.rmsd_ub:.3f}")
                if poses.crystal_rmsd is not None:
                    record.SetProp("crystal_rmsd", f"{poses.crystal_rmsd:.3f}")
            record.SetProp("docking_backend", result.backend)
            writer.write(record)
    finally:
        writer.close()
    if not target.is_file():
        raise DockFlowError(f"SDF writer produced no file: {target}")
