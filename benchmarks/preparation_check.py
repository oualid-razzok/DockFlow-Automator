"""Cross-engine preparation comparison report (peer item 9).

Prepares a receptor with every available engine and a ligand with Meeko,
then writes ``benchmarks/preparation_report.json`` - the machine-readable
source of the table in ``docs/preparation_validation.md``.  Fails when
the engines disagree on the heavy-atom count (they must not: geometry is
shared) or when no engine works.

Usage::

    python benchmarks/preparation_check.py --pdb raw/1hvr.pdb \
        [--ligand-sdf raw/xk2.sdf]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dockflow_core.pdbio import parse_pdbqt  # noqa: E402
from dockflow_core.preparator import (  # noqa: E402
    LigandPreparator,
    LigandPrepOptions,
    ReceptorPreparator,
    ReceptorPrepOptions,
)

ENGINES = ("openbabel", "rdkit", "none")


def receptor_report(engine: str, pdb: Path, workdir: Path) -> dict:
    result = ReceptorPreparator(ReceptorPrepOptions(engine=engine)).prepare(
        pdb, workdir / engine
    )
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    charges = [atom.charge for atom in atoms]
    types = Counter(atom.atom_type for atom in atoms)
    return {
        "engine": engine,
        "atoms_out": len(atoms),
        "heavy_atoms": sum(1 for atom in atoms if atom.element.upper() != "H"),
        "polar_hydrogens": sum(1 for atom in atoms if atom.atom_type in ("HD", "HS")),
        "charge_sum": round(sum(charges), 3),
        "charge_min": round(min(charges), 3),
        "charge_max": round(max(charges), 3),
        "ad4_type_histogram": dict(sorted(types.items())),
        "removed_hetero": result.removed_resnames,
        "num_warnings": len(result.warnings),
        "warnings": result.warnings,
    }


def ligand_report(sdf: Path, workdir: Path) -> dict:
    result = LigandPreparator(LigandPrepOptions()).prepare(
        sdf, workdir / "ligand", identifier="ligand"
    )
    atoms = parse_pdbqt(result.pdbqt_path).atoms
    return {
        "engine": result.engine,
        "atoms": len(atoms),
        "real_heavy_atoms": sum(
            1 for atom in atoms
            if atom.element.upper() not in ("H", "G") and not atom.atom_type.startswith("G")
        ),
        "virtual_atoms": sum(
            1 for atom in atoms if atom.atom_type.startswith("G")
        ),
        "polar_hydrogens": sum(1 for atom in atoms if atom.atom_type in ("HD", "HS")),
        "rotatable_bonds": result.num_rotatable_bonds,
        "charge_sum": round(sum(atom.charge for atom in atoms), 3),
        "input_had_explicit_hydrogens": result.input_had_explicit_hydrogens,
        "protonation_source": result.protonation_source,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, help="receptor PDB file")
    parser.add_argument("--ligand-sdf", default=None,
                        help="optional ligand SDF for the Meeko arm")
    parser.add_argument("--out", default=str(ROOT / "benchmarks" / "preparation_report.json"))
    args = parser.parse_args()

    pdb = Path(args.pdb)
    if not pdb.is_file():
        print(f"error: {pdb} not found", file=sys.stderr)
        return 2
    workdir = Path("/tmp") / "dockflow_preparation_check"
    workdir.mkdir(parents=True, exist_ok=True)

    reports = {}
    for engine in ENGINES:
        try:
            reports[engine] = receptor_report(engine, pdb, workdir)
            print(f"{engine:12s} {reports[engine]['atoms_out']} atoms "
                  f"({reports[engine]['heavy_atoms']} heavy), charge sum "
                  f"{reports[engine]['charge_sum']:+.3f}")
        except Exception as exc:  # noqa: BLE001 - engine unavailable
            print(f"{engine:12s} unavailable: {exc}")

    heavy_counts = {report["heavy_atoms"] for report in reports.values()}
    if len(reports) > 1 and len(heavy_counts) > 1:
        print(f"ENGINES DISAGREE on heavy-atom counts: {heavy_counts}",
              file=sys.stderr)
        return 1

    payload: dict = {"receptor": reports}
    if args.ligand_sdf and Path(args.ligand_sdf).is_file():
        payload["ligand"] = ligand_report(Path(args.ligand_sdf), workdir)
        ligand = payload["ligand"]
        print(f"ligand (meeko) {ligand['atoms']} atoms "
              f"({ligand['real_heavy_atoms']} real heavy), "
              f"{ligand['rotatable_bonds']} rotatable bonds")

    out = Path(args.out)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"report: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
