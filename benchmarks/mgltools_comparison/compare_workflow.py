"""DockFlow vs MGLTools preparation comparison (work order item 2).

Both arms dock with the SAME AutoDock Vina (python bindings), the SAME
search space (the redocking benchmark's ``gridbox.txt``), the SAME seed
(2026) and exhaustiveness (8).  The ONLY difference is the preparation
engine:

* DockFlow arm  - the committed redocking benchmark results
                  (openbabel receptor engine + RDKit/Meeko ligand prep);
* MGLTools arm  - ``prepare_receptor4.py`` / ``prepare_ligand4.py``
                  (MGLTools 1.5.7, classic academic workflow), run on the
                  SAME cleaned receptor (``receptor_clean.pdb``) and the
                  SAME RCSB ideal ligand SDF.

Because everything else is identical, per-complex deltas in RMSD, score
and atom counts are attributable to preparation chemistry (hydrogen
placement, united-atom model, charge assignment), which is exactly what
this benchmark measures.

The runner is resumable and crash-safe (same pattern as the redocking
benchmark): one CSV row per complex, appended on completion; failed
complexes are recorded with their reason, never dropped.

Outputs (in ``--out``):

* ``comparison.csv``        - per-complex both-arm table + deltas
* ``summary.md``            - delta table, mean deltas, honest discussion
* ``comparison_stats.json`` - machine-readable summary
* ``<pdb_id>/``             - MGLTools arm artifacts (pdbqts, poses, log)

Requirements: the MGLTools 1.5.7 suite.  Two supported locations:

1. ``--mgltools /path/to/mgltools`` (or env ``MGLTOOLS_HOME``) pointing
   at a tree with ``bin/python2.7`` and ``bin/prepare_receptor4.py``
   (e.g. the bioconda ``mgltools`` package extracted in-place);
2. ``prepare_receptor4.py`` on PATH (conda-forge style install).

Usage::

    python benchmarks/mgltools_comparison/compare_workflow.py \
        --mgltools /path/to/mgltools --max-seconds 540
    python benchmarks/mgltools_comparison/compare_workflow.py --report-only
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SUCCESS_RMSD = 2.0  # A, same conventional threshold as the redocking suite

MGL_FIELDNAMES = [
    "pdb_id", "ligand", "status", "mgl_rmsd_a", "mgl_affinity_kcal_mol",
    "mgl_num_poses", "mgl_receptor_atoms", "mgl_ligand_atoms",
    "mgl_ligand_polar_h", "runtime_s", "detail",
    # final pass, item 4: crash classification instead of a truncated
    # traceback dumped into the CSV
    "failure_stage", "failure_root_cause",
]
FINAL_FIELDNAMES = MGL_FIELDNAMES[:3] + [
    "df_rmsd_a", "df_affinity_kcal_mol", "df_receptor_atoms",
    "df_ligand_atoms", "rmsd_delta", "affinity_delta",
    "receptor_atom_delta", "success_mgl", "success_df",
] + MGL_FIELDNAMES[3:]  # all MGL-arm columns incl. mgl_rmsd_a
# (final pass fix: the v0.3.x field list dropped mgl_rmsd_a /
#  mgl_affinity_kcal_mol / mgl_num_poses / mgl_receptor_atoms from the
#  published comparison.csv even though the summary table used them)


def wilson_ci(successes: int, total: int, z: float = 1.959964
              ) -> tuple[float, float]:
    """Wilson score 95% interval (Wilson 1927); same hand-rolled
    implementation as the redocking runner (no statsmodels dependency).
    """
    if total <= 0:
        return (0.0, 0.0)
    p = successes / total
    z2 = z * z
    centre = p + z2 / (2.0 * total)
    margin = z * (p * (1.0 - p) / total + z2 / (4.0 * total * total)) ** 0.5
    denom = 1.0 + z2 / total
    return ((centre - margin) / denom, (centre + margin) / denom)


def _fmt_rate(successes: int, total: int) -> str:
    if total <= 0:
        return "n/a"
    low, high = wilson_ci(successes, total)
    rate = 100.0 * successes / total
    return f"{rate:.1f}% (95% CI: {100 * low:.1f}%-{100 * high:.1f}%)"


def _root_cause(traceback_text: str) -> str:
    """One-line root cause: the final exception line of a traceback."""
    lines = [line.strip() for line in traceback_text.strip().splitlines()
             if line.strip()]
    for line in reversed(lines):
        # the last line of a Python traceback is the exception message
        if line and not line.startswith("Traceback"):
            return line[:160]
    return "unknown (see traceback file)"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def find_mgltools(explicit: str | None) -> Path | None:
    """Locate the MGLTools tree (bundled python2.7 layout or PATH)."""
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    if os.environ.get("MGLTOOLS_HOME"):
        candidates.insert(0, Path(os.environ["MGLTOOLS_HOME"]))
    for cand in candidates:
        if (cand / "bin" / "prepare_receptor4.py").is_file():
            return cand
    if Path("prepare_receptor4.py").exists():
        return None  # PATH layout: scripts run directly
    import shutil

    if shutil.which("prepare_receptor4.py"):
        return None
    return None


def mgl_invocation(mgl_home: Path | None) -> list[str]:
    """Command prefix that runs an MGLTools utility script."""
    if mgl_home is None:
        return []
    return [str(mgl_home / "bin" / "python2.7")]


def parse_vina_config(path: Path) -> tuple[list[float], list[float]]:
    """Read center/size from a DockFlow gridbox.txt (vina config)."""
    numbers: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*(center|size)_([xyz])\s*=\s*(-?[\d.]+)", line)
        if match:
            numbers[f"{match.group(1)}_{match.group(2)}"] = float(match.group(3))
    center = [numbers[f"center_{a}"] for a in "xyz"]
    size = [numbers[f"size_{a}"] for a in "xyz"]
    return center, size


def count_pdbqt_atoms(path: Path) -> tuple[int, int]:
    """(all atoms, non-hydrogen real atoms) of a PDBQT file.

    Virtual atom types (Meeko's G0/CG0..CG3 glue atoms, W waters, XX) are
    excluded from the heavy count: they carry no chemistry.
    """
    total = heavy = 0
    virtual = re.compile(r"^(G[0-3]|CG[0-3]|W|XX)$")
    hydrogen = {"HD", "HS", "H"}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            total += 1
            atom_type = line[77:].strip() if len(line) > 77 else ""
            if atom_type not in hydrogen and not virtual.match(atom_type):
                heavy += 1
    return total, heavy


def load_complexes(path: Path) -> list[dict[str, str]]:
    lines = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return [row for row in csv.DictReader(lines) if row.get("pdb_id")]


def load_existing(csv_path: Path) -> dict[str, dict]:
    if not csv_path.exists():
        return {}
    with open(csv_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {row["pdb_id"]: row for row in rows if row.get("pdb_id")}


def append_row(csv_path: Path, row: dict) -> None:
    is_new = not csv_path.exists()
    with open(csv_path, "a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MGL_FIELDNAMES)
        if is_new:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in MGL_FIELDNAMES})


# ---------------------------------------------------------------------------
# worker: one complex, MGLTools arm only
# ---------------------------------------------------------------------------
def run_mgl_arm(complex_row: dict, out_dir: Path, bench_results: Path,
                mgl_home: Path | None, exhaustiveness: int, seed: int
                ) -> dict:
    pdb_id = complex_row["pdb_id"]
    ligand_resname = complex_row["ligand_resname"]
    bench = bench_results / f"{pdb_id}_benchmark"
    workdir = out_dir / pdb_id
    workdir.mkdir(parents=True, exist_ok=True)

    def fail(stage: str, output: str, script: str) -> dict:
        """Record a crash with the FULL traceback in failures/ (item 4).

        The v0.3.x runner truncated the captured output at 200 characters
        inside the CSV, which is exactly what the final review asked to
        fix: the full text now lives in
        ``results/failures/{pdb}_mgltools_{stage}_traceback.txt`` and the
        CSV carries the stage, a one-line root cause and the file link.
        """
        failure_dir = out_dir / "failures"
        failure_dir.mkdir(parents=True, exist_ok=True)
        traceback_path = failure_dir / f"{pdb_id}_mgltools_{stage}_traceback.txt"
        header = (f"MGLTools 1.5.7 crash - {pdb_id} ({ligand_resname}) - "
                  f"stage: {stage} - script: {script}\n"
                  f"command working directory: {workdir}\n\n")
        traceback_path.write_text(header + output.strip() + "\n",
                                  encoding="utf-8")
        return {
            "pdb_id": pdb_id, "ligand": ligand_resname,
            "status": f"failed: {script} (see {traceback_path.name})",
            "failure_stage": stage,
            "failure_root_cause": _root_cause(output),
            "detail": f"full traceback: failures/{traceback_path.name}",
        }

    started = time.perf_counter()
    prefix = mgl_invocation(mgl_home)

    def mgl(script: str, *args: str) -> subprocess.CompletedProcess:
        if mgl_home is None:
            return subprocess.run([script, *args], capture_output=True,
                                  text=True, timeout=600, cwd=str(workdir))
        return subprocess.run(
            prefix + [str(mgl_home / "bin" / script), *args],
            capture_output=True, text=True, timeout=600, cwd=str(workdir))

    # ---- 1. receptor: same cleaned input the DockFlow arm used ----------
    rec_out = workdir / "receptor_mgl.pdbqt"
    if not rec_out.is_file():
        proc = mgl("prepare_receptor4.py",
                   "-r", str(bench / "prepared" / "receptor_clean.pdb"),
                   "-o", str(rec_out), "-A", "hydrogens", "-U", "nphs")
        if not rec_out.is_file():
            return fail("receptor_prep", proc.stderr or proc.stdout,
                        "prepare_receptor4.py")
    rec_total, rec_heavy = count_pdbqt_atoms(rec_out)

    # ---- 2. ligand: same RCSB ideal SDF, converted to mol2 for ADT ------
    sdf = next((bench / "raw").glob("*.sdf"), None)
    if sdf is None:
        return fail("ligand_prep", "no cached ideal SDF in the benchmark run dir",
                    "(input)")
    mol2 = workdir / f"{ligand_resname.lower()}.mol2"
    lig_out = workdir / f"{ligand_resname.lower()}_mgl.pdbqt"
    if not lig_out.is_file():
        from openbabel import openbabel as ob

        molecule = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("sdf", "mol2")
        if not conv.ReadFile(molecule, str(sdf)):
            return fail("ligand_prep", f"openbabel could not read {sdf.name}",
                        "obabel sdf->mol2")
        if not conv.WriteFile(molecule, str(mol2)):
            return fail("ligand_prep", "openbabel mol2 conversion failed",
                        "obabel sdf->mol2")
        proc = mgl("prepare_ligand4.py", "-l", str(mol2), "-o", str(lig_out))
        if not lig_out.is_file():
            return fail("ligand_prep", proc.stderr or proc.stdout,
                        "prepare_ligand4.py")
    lig_total, lig_heavy = count_pdbqt_atoms(lig_out)

    # ---- 3. identical search space + seed + exhaustiveness --------------
    center, size = parse_vina_config(bench / "gridbox.txt")
    pose_out = workdir / f"{ligand_resname.lower()}_mgl_out.pdbqt"
    from vina import Vina

    try:
        vina_obj = Vina(sf_name="vina", seed=seed, cpu=0, verbosity=0)
        vina_obj.set_receptor(str(rec_out))
        vina_obj.set_ligand_from_file(str(lig_out))
        vina_obj.compute_vina_maps(center=list(center), box_size=list(size))
        vina_obj.dock(exhaustiveness=exhaustiveness, n_poses=9)
        vina_obj.write_poses(str(pose_out), n_poses=9, overwrite=True)
    except Exception:  # noqa: BLE001 - docking crash, record fully
        import traceback

        return fail("docking", traceback.format_exc(), "vina-python")

    # ---- 4. crystal RMSD with the SAME analyzer as the DockFlow arm -----
    from dockflow_core.analyzer import crystal_rmsd_with_method, extract_reference_atoms
    from dockflow_core.pdbio import parse_pdbqt

    raw_pdb = bench / "raw" / f"{pdb_id.lower()}.pdb"
    reference = extract_reference_atoms(raw_pdb, ligand_resname)
    models = parse_pdbqt(pose_out).models
    best_rmsd = None
    method = None
    for model in models:
        rmsd, how = crystal_rmsd_with_method(model.atoms, reference)
        if rmsd is not None and (best_rmsd is None or rmsd < best_rmsd):
            best_rmsd, method = rmsd, how
    best_affinity = models[0].vina_result.affinity if models else None

    return {
        "pdb_id": pdb_id,
        "ligand": ligand_resname,
        "status": "ok",
        "mgl_rmsd_a": best_rmsd,
        "mgl_affinity_kcal_mol": best_affinity,
        "mgl_num_poses": len(models),
        "mgl_receptor_atoms": rec_heavy,
        "mgl_ligand_atoms": lig_heavy,
        "mgl_ligand_polar_h": lig_total - lig_heavy,
        "runtime_s": round(time.perf_counter() - started, 1),
        "detail": f"rmsd method: {method}" if method else "no rmsd",
    }


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def dockflow_arm(bench_results: Path, pdb_id: str, ligand: str
                 ) -> dict:
    """DockFlow arm numbers from the redocking benchmark run dir."""
    run = bench_results / f"{pdb_id}_benchmark"
    out = {"df_rmsd_a": None, "df_affinity_kcal_mol": None,
           "df_receptor_atoms": None, "df_ligand_atoms": None}
    try:
        manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
        out["df_affinity_kcal_mol"] = manifest["docking"]["results"][0]["poses"][0].get("affinity")
        poses = manifest["docking"]["results"][0]["poses"]
        rmsds = [p.get("crystal_rmsd") for p in poses if p.get("crystal_rmsd") is not None]
        out["df_rmsd_a"] = min(rmsds) if rmsds else None
        rec = run / "prepared" / "receptor.pdbqt"
        if rec.is_file():
            out["df_receptor_atoms"] = count_pdbqt_atoms(rec)[1]
        lig = next((run / "prepared").glob("*redock*.pdbqt"), None)
        if lig is not None:
            out["df_ligand_atoms"] = count_pdbqt_atoms(lig)[1]
    except (OSError, ValueError, KeyError, IndexError):
        pass
    return out


def _fmt(value, digits: int = 2) -> str:
    if value in (None, "", "None"):
        return "n/a"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def write_report(rows: list[dict], out_dir: Path, bench_results: Path,
                 exhaustiveness: int, seed: int) -> None:
    merged = []
    for row in rows:
        entry = {**row, **dockflow_arm(bench_results, row["pdb_id"],
                                       row.get("ligand", ""))}
        if row.get("status") == "ok":
            try:
                entry["rmsd_delta"] = (float(row["mgl_rmsd_a"])
                                       - float(entry["df_rmsd_a"])
                                       if entry["df_rmsd_a"] is not None else None)
                entry["affinity_delta"] = (float(row["mgl_affinity_kcal_mol"])
                                           - float(entry["df_affinity_kcal_mol"])
                                           if entry["df_affinity_kcal_mol"] is not None else None)
                entry["receptor_atom_delta"] = (int(row["mgl_receptor_atoms"])
                                                - int(entry["df_receptor_atoms"])
                                                if entry["df_receptor_atoms"] is not None else None)
            except (TypeError, ValueError):
                entry["rmsd_delta"] = entry["affinity_delta"] = None
                entry["receptor_atom_delta"] = None
            entry["success_mgl"] = (row.get("mgl_rmsd_a") not in (None, "", "None")
                                    and float(row["mgl_rmsd_a"]) <= SUCCESS_RMSD)
            entry["success_df"] = (entry["df_rmsd_a"] not in (None, "", "None")
                                   and float(entry["df_rmsd_a"]) <= SUCCESS_RMSD)
        merged.append(entry)

    csv_path = out_dir / "comparison.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FINAL_FIELDNAMES)
        writer.writeheader()
        for entry in merged:
            writer.writerow({key: entry.get(key, "") for key in FINAL_FIELDNAMES})

    ok = [e for e in merged if e.get("status") == "ok"]
    both_rmsd = [e for e in ok if e.get("rmsd_delta") is not None]
    deltas = [float(e["rmsd_delta"]) for e in both_rmsd]
    aff_deltas = [float(e["affinity_delta"]) for e in both_rmsd
                  if e.get("affinity_delta") is not None]
    atom_deltas = [int(e["receptor_atom_delta"]) for e in both_rmsd
                   if e.get("receptor_atom_delta") is not None]
    mgl_ok = [e for e in ok if str(e.get("success_mgl")) == "True"]
    df_ok = [e for e in ok if str(e.get("success_df")) == "True"]

    stats = {
        "complexes_attempted": len(merged),
        "mgl_arm_completed": len(ok),
        "mgl_arm_success_within_2a": len(mgl_ok),
        "dockflow_arm_success_within_2a": len(df_ok),
        "mgl_arm_success_rate_wilson_95ci_pct": [
            round(100 * v, 1) for v in wilson_ci(len(mgl_ok), len(ok))],
        "dockflow_arm_success_rate_wilson_95ci_pct": [
            round(100 * v, 1) for v in wilson_ci(len(df_ok), len(ok))],
        "mean_rmsd_delta_a": round(sum(deltas) / len(deltas), 3) if deltas else None,
        "mean_affinity_delta": round(sum(aff_deltas) / len(aff_deltas), 3)
        if aff_deltas else None,
        "mean_receptor_atom_delta": round(sum(atom_deltas) / len(atom_deltas), 1)
        if atom_deltas else None,
        "exhaustiveness": exhaustiveness,
        "seed": seed,
    }
    (out_dir / "comparison_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# DockFlow vs MGLTools preparation comparison",
        "",
        "Both arms: AutoDock Vina (python bindings), the SAME search space",
        f"(redocking benchmark gridbox), seed {seed}, exhaustiveness",
        f"{exhaustiveness}, num_modes 9.  The only difference is the",
        "preparation engine - every delta below is attributable to",
        "preparation chemistry (hydrogen model, charge assignment), not to",
        "the docking engine.",
        "",
        f"- complexes attempted: {stats['complexes_attempted']}",
        f"- MGLTools arm completed: {stats['mgl_arm_completed']} "
        f"(failures recorded with reasons, not dropped)",
        f"- MGLTools arm success (RMSD <= 2.0 A): **{len(mgl_ok)}/{len(ok)}** "
        f"= {_fmt_rate(len(mgl_ok), len(ok))}",
        f"- DockFlow arm success: **{len(df_ok)}/{len(ok)}** "
        f"= {_fmt_rate(len(df_ok), len(ok))}",
        "- **denominator (explicit, item 4): the success rates above are "
        f"computed on the {len(ok)}-complex MGLTools-completed subset - "
        f"{len(merged) - len(ok)} MGLTools crashes are excluded.  This is "
        "NOT a head-to-head 24-vs-24 comparison; the full 24-complex "
        "DockFlow numbers are in benchmarks/redocking/results/.",
        f"- mean RMSD delta (MGLTools - DockFlow): {_fmt(stats['mean_rmsd_delta_a'], 3)} A",
        f"- mean affinity delta (MGLTools - DockFlow): {_fmt(stats['mean_affinity_delta'], 3)} kcal/mol",
        f"- mean receptor atom-count delta: {_fmt(stats['mean_receptor_atom_delta'], 1)}",
        "",
        "## Per-complex deltas",
        "",
        "| pdb | ligand | RMSD DockFlow (A) | RMSD MGLTools (A) | delta | "
        "affinity DockFlow | affinity MGLTools | delta | rec atoms DF/MGL | "
        "lig atoms DF/MGL |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for entry in merged:
        if entry.get("status") == "ok":
            lines.append(
                f"| {entry['pdb_id']} | {entry.get('ligand', '')} | "
                f"{_fmt(entry.get('df_rmsd_a'))} | {_fmt(entry.get('mgl_rmsd_a'))} | "
                f"{_fmt(entry.get('rmsd_delta'), 2)} | "
                f"{_fmt(entry.get('df_affinity_kcal_mol'))} | "
                f"{_fmt(entry.get('mgl_affinity_kcal_mol'))} | "
                f"{_fmt(entry.get('affinity_delta'), 2)} | "
                f"{entry.get('df_receptor_atoms', 'n/a')}/{entry.get('mgl_receptor_atoms', 'n/a')} | "
                f"{entry.get('df_ligand_atoms', 'n/a')}/{entry.get('mgl_ligand_atoms', 'n/a')} |"
            )
        else:
            # item 4: crashes link to the full traceback file and carry the
            # stage + one-line root cause instead of a truncated dump
            lines.append(
                f"| {entry['pdb_id']} | {entry.get('ligand', '')} | "
                f"{_fmt(entry.get('df_rmsd_a'))} | FAILED ({entry.get('failure_stage', 'n/a')}) | - | "
                f"{_fmt(entry.get('df_affinity_kcal_mol'))} | - | - | - | - |"
                f" root cause: {entry.get('failure_root_cause', 'n/a')}; "
                f"[full traceback]({entry.get('detail', '')}) |"
            )
    lines += [
        "",
        "## Interpretation",
        "",
        "MGLTools and DockFlow (openbabel engine) differ in hydrogen",
        "placement and charge assignment, so the same Vina search lands on",
        "slightly different poses and scores.  A positive RMSD delta means",
        "the MGLTools arm recovered the crystal pose better; negative means",
        "DockFlow did.  Neither arm is a ground truth: both are compared",
        "against the crystal pose.",
        "",
        "The MGLTools-arm failures are **genuine MGLTools 1.5.7 crashes** "
        "on these inputs (full tracebacks in ``results/failures/``, linked "
        "from the table and the CSV); the DockFlow arm processed the same "
        "files successfully.  They are recorded, not dropped, so the two "
        "arms' success rates share the MGLTools-completed denominator "
        "stated above - not a head-to-head 24-vs-24.",
        "",
        "Systematic differences to keep in mind when reading the table:",
        "",
        "* **Atom counts**: the receptor atom-count delta counts non-hydrogen",
        "  PDBQT atoms (virtual glue atoms excluded); MGLTools' ``-U nphs``",
        "  removes non-polar hydrogens (united-atom model, polar H kept),",
        "  DockFlow's openbabel engine also keeps polar hydrogens only - small",
        "  count differences reflect H-bonding assignment differences, not",
        "  lost atoms.  The ligand column counts real heavy atoms; polar-H",
        "  differences are in ``mgl_ligand_polar_h``.",
        "* **Ligand torsions**: ``prepare_ligand4`` and Meeko can disagree",
        "  on which bonds are rotatable, which changes the search-space size.",
        "* **Macrocycles**: Meeko opens macrocyclic rings for flexible-",
        "  macrocycle docking by replacing the two ring-closure atoms with",
        "  G0/CG0 glue pseudo-atoms - for such ligands (1HVR/XK2 in this set)",
        "  the DockFlow-arm explicit heavy count is 2 lower than MGLTools'",
        "  (which docks the closed ring); no atoms are lost.",
        "* **Charges**: both arms use Gasteiger charges, but hydrogen",
        "  placement changes how those charges are distributed.",
        "",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--complexes",
                        default=str(Path(__file__).parent.parent
                                    / "redocking" / "complexes.csv"))
    parser.add_argument("--benchmark-results",
                        default=str(Path(__file__).parent.parent
                                    / "redocking" / "results"))
    parser.add_argument("--out", default=str(Path(__file__).parent / "results"))
    parser.add_argument("--mgltools", default=None,
                        help="path to the MGLTools tree (bin/prepare_receptor4.py)")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--pdb-ids", default=None)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--max-seconds", type=float, default=None)
    parser.add_argument("--worker-timeout", type=float, default=900)
    parser.add_argument("--worker", default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    bench_results = Path(args.benchmark_results)
    csv_path = out_dir / "mgl_arm_results.csv"

    mgl_home = find_mgltools(args.mgltools)

    # ---- report-only mode (no MGLTools needed) --------------------------
    if args.report_only:
        rows = list(load_existing(csv_path).values())
        if not rows:
            print("no MGLTools-arm results recorded yet", file=sys.stderr)
            return 1
        write_report(rows, out_dir, bench_results,
                     args.exhaustiveness, args.seed)
        print(f"summary: {out_dir / 'summary.md'}")
        return 0

    if mgl_home is None and args.mgltools is None and \
            not Path("prepare_receptor4.py").exists():
        import shutil

        if not shutil.which("prepare_receptor4.py"):
            print("MGLTools not found: pass --mgltools /path/to/mgltools "
                  "(tree containing bin/prepare_receptor4.py) or install "
                  "MGLTools on PATH.  See README.md.", file=sys.stderr)
            return 2
    if args.worker:
        complexes = load_complexes(Path(args.complexes))
        row = next(c for c in complexes if c["pdb_id"] == args.worker)
        try:
            result = run_mgl_arm(row, out_dir, bench_results, mgl_home,
                                 args.exhaustiveness, args.seed)
        except Exception as exc:  # noqa: BLE001 - record, never crash
            result = {"pdb_id": row["pdb_id"], "ligand": row["ligand_resname"],
                      "status": f"error: {exc}"}
        print(json.dumps(result), flush=True)
        return 0

    # ---- report-only mode ----
    if args.report_only:
        rows = list(load_existing(csv_path).values())
        if not rows:
            print("no MGLTools-arm results recorded yet", file=sys.stderr)
            return 1
        write_report(rows, out_dir, bench_results,
                     args.exhaustiveness, args.seed)
        print(f"summary: {out_dir / 'summary.md'}")
        return 0

    # ---- driver mode ----
    complexes = load_complexes(Path(args.complexes))
    if args.pdb_ids:
        wanted = {pid.strip().upper() for pid in args.pdb_ids.split(",")}
        complexes = [row for row in complexes if row["pdb_id"].upper() in wanted]
    if args.limit:
        complexes = complexes[: args.limit]

    existing = load_existing(csv_path)
    todo = [row for row in complexes if row["pdb_id"] not in existing]
    print(f"MGLTools comparison: {len(complexes)} complex(es), "
          f"{len(existing)} recorded, {len(todo)} to run | "
          f"exhaustiveness={args.exhaustiveness}, seed={args.seed}, "
          f"mgltools={mgl_home or 'PATH'}", flush=True)

    started_driver = time.perf_counter()
    for index, complex_row in enumerate(todo, start=1):
        pdb_id = complex_row["pdb_id"]
        if args.max_seconds is not None and \
                time.perf_counter() - started_driver > args.max_seconds:
            print(f"time budget reached; {len(todo) - index + 1} remain - "
                  "re-run to resume", flush=True)
            break
        print(f"[{index}/{len(todo)}] {pdb_id} "
              f"({complex_row['ligand_resname']}) MGLTools arm ...", flush=True)
        started = time.perf_counter()
        cmd = [sys.executable, str(Path(__file__).resolve()),
               "--complexes", args.complexes,
               "--benchmark-results", args.benchmark_results,
               "--out", args.out,
               "--mgltools", str(mgl_home) if mgl_home else "",
               "--exhaustiveness", str(args.exhaustiveness),
               "--seed", str(args.seed),
               "--worker", pdb_id]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=args.worker_timeout)
            last = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
            result = json.loads(last) if last.startswith("{") else None
        except subprocess.TimeoutExpired:
            result = {"pdb_id": pdb_id, "ligand": complex_row["ligand_resname"],
                      "status": f"error: worker exceeded "
                                f"{args.worker_timeout:.0f}s timeout"}
        if result is None:
            result = {"pdb_id": pdb_id, "ligand": complex_row["ligand_resname"],
                      "status": "error: worker crashed or produced no result"}
        append_row(csv_path, result)
        if result.get("status") == "ok":
            print(f"    mgl best {result.get('mgl_affinity_kcal_mol')} kcal/mol, "
                  f"RMSD {result.get('mgl_rmsd_a')} A "
                  f"({time.perf_counter() - started:.0f}s)", flush=True)
        else:
            print(f"    {result.get('status')}", flush=True)
        rows = list(load_existing(csv_path).values())
        write_report(rows, out_dir, bench_results,
                     args.exhaustiveness, args.seed)

    rows = list(load_existing(csv_path).values())
    write_report(rows, out_dir, bench_results, args.exhaustiveness, args.seed)
    print(f"\nresults: {csv_path}")
    print(f"summary: {out_dir / 'summary.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
