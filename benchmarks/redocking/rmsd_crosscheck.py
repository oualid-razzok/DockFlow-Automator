"""Independent cross-validation of the symmetry-aware crystal RMSD (item 5).

DockFlow's published crystal RMSDs are computed by
``dockflow_core.analyzer.crystal_rmsd_with_method`` (RDKit ``GetBestRMS``
minimised over molecular-graph automorphisms).  This script validates
that implementation against a SECOND, independent implementation on
every pose of the 24-complex redocking benchmark:

* **primary**: `spyrmsd <https://github.com/RMeli/spyrmsd>`_ (BSD-2) - a
  different code base by different authors, called through its
  ``symmrmsd`` API (networkx graph-automorphism enumeration + its own
  QCP/Kabsch minimisation).  spyrmsd is a BENCHMARK-ONLY tool here, not
  a package dependency;
* **fallback** (when spyrmsd is not installed, e.g. in CI): a
  hand-rolled Kabsch + automorphism-enumeration implementation -
  automorphisms enumerated with RDKit ``GetSubstructMatches`` on a
  graph built by THIS script, RMSD minimised with THIS script's own
  numpy Kabsch over every automorphism permutation (the same algorithm
  spyrmsd implements, independently written).

The method actually used is printed, written to
``crosscheck_stats.json`` and recorded on the plot, so the check is
never silently method-less.

Inputs: the committed benchmark artifacts (per-pose PDBQT files +
co-crystallised reference ligands).  The atom-inclusion rule (heavy,
non-virtual atoms - meeko G0/CG0 glue atoms excluded) matches the
production analyzer so the two implementations receive identical atom
sets; everything after that (graph perception, automorphism search,
superposition) is independent.

Outputs (in ``--out``, default: results/):

* ``rmsd_crosscheck.png``    - scatter of DockFlow RMSD vs independent
  RMSD per pose, with the y=x line
* ``crosscheck_stats.json``  - per-pose deltas, max |delta|, method
* ``crosscheck.csv``         - one row per pose (pdb, pose, both RMSDs)

Exit code: non-zero when max |delta| >= 0.1 A (the item-5 acceptance
criterion) or when fewer than 100 poses could be compared.

Usage::

    python benchmarks/redocking/rmsd_crosscheck.py
    python benchmarks/redocking/rmsd_crosscheck.py --pdb-ids 1HVR,1STP
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dockflow_core.analyzer import (  # noqa: E402
    chemical_element,
    crystal_rmsd_with_method,
    extract_reference_atoms,
)
from dockflow_core.pdbio import parse_pdbqt  # noqa: E402

DEFAULT_RESULTS = Path(__file__).parent / "results"
DELTA_LIMIT_A = 0.1  # item 5: "assert max delta < 0.1 A"
MIN_POSES = 100

# Covalent radii (A) for the elements the benchmark ligands contain; a
# bond exists when d(i, j) < r_i + r_j + 0.40 A.  Used to build the
# adjacency matrices for spyrmsd and the fallback graph.
_COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "SI": 1.11, "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39,
}
_ATOMIC_NUMBERS = {
    "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "SI": 14, "P": 15, "S": 16,
    "CL": 17, "BR": 35, "I": 53,
}


def _heavy_atoms(atoms) -> list:
    """Production atom-inclusion rule: heavy, non-virtual atoms only."""
    return [a for a in atoms if chemical_element(a) not in (None, "H")]


def _coords(atoms) -> np.ndarray:
    return np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)


def _aprops(atoms) -> np.ndarray:
    return np.array([
        _ATOMIC_NUMBERS[chemical_element(a).upper()] for a in atoms
    ])


def _adjacency(coords: np.ndarray, aprops: np.ndarray) -> np.ndarray:
    n = len(coords)
    radii = np.array([
        _COVALENT_RADII[{v: k for k, v in _ATOMIC_NUMBERS.items()}[a]]
        for a in aprops
    ])
    am = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(coords[i] - coords[j]) < radii[i] + radii[j] + 0.40:
                am[i, j] = am[j, i] = 1
    return am


# ---------------------------------------------------------------------------
# Independent implementations
# ---------------------------------------------------------------------------
def _spyrmsd_rmsd(ref_atoms, pose_atoms):
    """spyrmsd symmrmsd (primary independent implementation)."""
    from spyrmsd.rmsd import symmrmsd

    coords_ref = _coords(ref_atoms)
    coords_pose = _coords(pose_atoms)
    aprops_ref = _aprops(ref_atoms)
    aprops_pose = _aprops(pose_atoms)
    am_ref = _adjacency(coords_ref, aprops_ref)
    am_pose = _adjacency(coords_pose, aprops_pose)
    return float(symmrmsd(coords_ref, coords_pose, aprops_ref, aprops_pose,
                          am_ref, am_pose, minimize=True))


def _handrolled_rmsd(ref_atoms, pose_atoms):
    """Hand-rolled Kabsch + automorphism enumeration (fallback).

    Automorphisms are enumerated with RDKit substructure self-matches on
    a molecule built by THIS function; the Kabsch superposition and the
    minimum over permutations are implemented here with numpy.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign  # noqa: F401 - documented contrast

    def _mol(atoms):
        n = len(atoms)
        mol = Chem.RWMol()
        conf = Chem.Conformer(n)
        for i, atom in enumerate(atoms):
            a = Chem.Atom(chemical_element(atom))
            mol.AddAtom(a)
            conf.SetAtomPosition(i, (atom.x, atom.y, atom.z))
        coords = _coords(atoms)
        aprops = _aprops(atoms)
        am = _adjacency(coords, aprops)
        for i in range(n):
            for j in range(i + 1, n):
                if am[i, j]:
                    mol.AddBond(i, j, Chem.BondType.SINGLE)
        mol = mol.GetMol()
        mol.AddConformer(conf, assignId=True)
        Chem.SanitizeMol(mol, catchErrors=True)
        return mol

    ref_mol = _mol(ref_atoms)
    automorphisms = ref_mol.GetSubstructMatches(
        ref_mol, uniquify=False, useChirality=False, maxMatches=10000)
    coords_ref = _coords(ref_atoms)
    coords_pose = _coords(pose_atoms)

    def _kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
        pc = p.mean(axis=0)
        qc = q.mean(axis=0)
        covariance = (p - pc).T @ (q - qc)
        u, _s, vt = np.linalg.svd(covariance)
        d = np.sign(np.linalg.det(u @ vt))
        rotation = u @ np.diag([1.0, 1.0, d]) @ vt
        aligned = (p - pc) @ rotation + qc
        return float(np.sqrt(np.mean(np.sum((aligned - q) ** 2, axis=1))))

    best = np.inf
    for automorphism in automorphisms:
        permuted = coords_pose[list(automorphism)]
        best = min(best, _kabsch_rmsd(permuted, coords_ref))
    return float(best)


def _independent_rmsd(ref_atoms, pose_atoms):
    """(rmsd, method name) via spyrmsd when importable, else hand-rolled."""
    try:
        return _spyrmsd_rmsd(ref_atoms, pose_atoms), "spyrmsd-symmrmsd"
    except ImportError:
        return _handrolled_rmsd(ref_atoms, pose_atoms), \
            "handrolled-kabsch-automorphism"


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def load_complexes(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    return [row for row in csv.DictReader(lines) if row.get("pdb_id")]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default=str(DEFAULT_RESULTS))
    parser.add_argument("--complexes",
                        default=str(Path(__file__).parent / "complexes.csv"))
    parser.add_argument("--out", default=str(DEFAULT_RESULTS))
    parser.add_argument("--pdb-ids", default=None)
    args = parser.parse_args()

    results_dir = Path(args.results)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    complexes = load_complexes(Path(args.complexes))
    if args.pdb_ids:
        wanted = {pid.strip().upper() for pid in args.pdb_ids.split(",")}
        complexes = [row for row in complexes
                     if row["pdb_id"].upper() in wanted]

    records: list[dict] = []
    method_used = None
    for complex_row in complexes:
        pdb_id = complex_row["pdb_id"]
        resname = complex_row["ligand_resname"]
        run_dir = results_dir / f"{pdb_id}_benchmark"
        raw_pdb = run_dir / "raw" / f"{pdb_id.lower()}.pdb"
        pose_files = sorted(run_dir.glob("docking/*_redock_out.pdbqt"))
        if not raw_pdb.is_file() or not pose_files:
            print(f"{pdb_id}: no committed poses/reference - skipped",
                  file=sys.stderr)
            continue
        try:
            reference = _heavy_atoms(extract_reference_atoms(raw_pdb, resname))
        except ValueError as exc:
            print(f"{pdb_id}: reference unavailable ({exc})", file=sys.stderr)
            continue
        for pose_file in pose_files:
            for model in parse_pdbqt(pose_file).models:
                pose_atoms = _heavy_atoms(model.atoms)
                dockflow, _method = crystal_rmsd_with_method(
                    pose_atoms, reference)
                if dockflow is None:
                    continue
                independent, method_used = _independent_rmsd(
                    reference, pose_atoms)
                records.append({
                    "pdb_id": pdb_id,
                    "pose": model.index,
                    "dockflow_rmsd_a": round(dockflow, 6),
                    "independent_rmsd_a": round(independent, 6),
                    "delta_a": round(independent - dockflow, 6),
                    "independent_method": method_used,
                })

    if not records:
        print("no poses could be cross-checked", file=sys.stderr)
        return 1
    deltas = [abs(r["delta_a"]) for r in records]
    max_delta = max(deltas)
    stats = {
        "poses_compared": len(records),
        "complexes": len({r["pdb_id"] for r in records}),
        "independent_method": method_used,
        "max_abs_delta_a": round(max_delta, 6),
        "mean_abs_delta_a": round(sum(deltas) / len(deltas), 6),
        "delta_limit_a": DELTA_LIMIT_A,
        "passed": max_delta < DELTA_LIMIT_A,
    }
    (out_dir / "crosscheck_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    with open(out_dir / "crosscheck.csv", "w", newline="",
              encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    # scatter: DockFlow vs independent, with y=x line
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        xs = [r["independent_rmsd_a"] for r in records]
        ys = [r["dockflow_rmsd_a"] for r in records]
        fig, ax = plt.subplots(figsize=(6.0, 6.0), constrained_layout=True)
        ax.scatter(xs, ys, s=18, alpha=0.65, color="#4c72b0", zorder=3)
        lim = max(max(xs), max(ys)) * 1.05
        ax.plot([0, lim], [0, lim], color="#c62828", ls="--", lw=1.2,
                label="y = x (perfect agreement)", zorder=2)
        ax.set_xlabel(f"independent RMSD ({method_used}) (A)")
        ax.set_ylabel("DockFlow RMSD (rdkit-getbestrms) (A)")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_title(f"Symmetry-aware RMSD cross-check: {len(records)} poses, "
                     f"max |delta| = {max_delta:.4f} A")
        ax.legend(loc="lower right")
        fig.savefig(out_dir / "rmsd_crosscheck.png", dpi=150)
        plt.close(fig)
    except ImportError:
        print("matplotlib unavailable - scatter skipped", file=sys.stderr)

    print(json.dumps(stats, indent=2))
    if not stats["passed"]:
        print(f"FAIL: max |delta| {max_delta:.4f} A >= "
              f"{DELTA_LIMIT_A} A", file=sys.stderr)
        return 1
    if stats["poses_compared"] < MIN_POSES:
        print(f"FAIL: only {stats['poses_compared']} poses compared "
              f"(< {MIN_POSES})", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
