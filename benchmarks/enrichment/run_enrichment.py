"""DUD-E enrichment benchmark runner (work order item 3).

Docks a DUD-E actives + decoys panel against the DockFlow example
receptor (HIV-1 protease, 1HVR) with the ligand-derived grid box and
publishes the enrichment metrics (ROC AUC, EF@1%, EF@5%, BEDROC).

The runner is resumable like the redocking benchmark:

* one CSV row per ligand, appended the moment it is scored,
* prepared ligand PDBQTs are cached (re-runs skip them),
* ``--report-only`` recomputes metrics/plots from the CSV,
* ``--max-seconds`` stops cleanly between ligands for chunked execution.

Protocol (documented honestly in the results):

* target: DUD-E ``hivpr`` (consistent with the 1HVR example),
* panel: seeded subsample of the DUD-E actives/decoys ISM lists
  (the full 536 actives / 35,750 decoys set is far beyond a single-CPU
  benchmark; the subsample ratio is recorded in the summary),
* exhaustiveness 4, num_modes 1 (screening-grade settings: enrichment
  depends on the BEST pose score only),
* receptor + box: the 1HVR redocking benchmark artifacts (box derived
  from the XK2 co-crystal ligand - strong assumption, documented).

Usage::

    python benchmarks/enrichment/run_enrichment.py --max-seconds 540
    python benchmarks/enrichment/run_enrichment.py --report-only
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

DEFAULT_DATA = Path(__file__).parent / "data"
DEFAULT_REF = Path(__file__).resolve().parents[1] / "redocking" / "results" / "1HVR_benchmark"

FIELDNAMES = ["ligand", "label", "status", "best_affinity_kcal_mol", "runtime_s"]


def load_panel(data_dir: Path, n_actives: int, n_decoys: int,
               seed: int) -> list[dict]:
    """Seeded subsample of the DUD-E ISM lists: [{id, smiles, label}]."""
    rng = random.Random(seed)

    def parse(path: Path, label: bool) -> list[dict]:
        rows = []
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 2:
                rows.append({"id": parts[1], "smiles": parts[0],
                             "label": label})
        return rows

    actives = parse(data_dir / "actives_final.ism", True)
    decoys = parse(data_dir / "decoys_final.ism", False)
    if n_actives < len(actives):
        actives = rng.sample(actives, n_actives)
    if n_decoys < len(decoys):
        decoys = rng.sample(decoys, n_decoys)
    return actives + decoys


def load_existing(csv_path: Path) -> dict[str, dict]:
    if not csv_path.exists():
        return {}
    with open(csv_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {row["ligand"]: row for row in rows if row.get("ligand")}


def append_row(csv_path: Path, row: dict) -> None:
    is_new = not csv_path.exists()
    with open(csv_path, "a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if is_new:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in FIELDNAMES})


def run_one(entry: dict, out_dir: Path, receptor: Path, box: tuple,
            exhaustiveness: int, seed: int) -> dict:
    """Prepare (cached) + dock one panel ligand; return its CSV row."""
    from dockflow_core.docker_engine import VinaConfig, VinaEngine
    from dockflow_core.preparator import LigandPreparator, LigandPrepOptions

    ligand_id = entry["id"]
    pdbqt = out_dir / "prepared" / f"{ligand_id}.pdbqt"
    started = time.perf_counter()
    status = "ok"
    affinity = None
    if not pdbqt.is_file():
        preparator = LigandPreparator(LigandPrepOptions(
            embed_3d=True, minimize=True, random_seed=seed))
        try:
            result = preparator.prepare(entry["smiles"],
                                        out_dir / "prepared", ligand_id)
            pdbqt = result.pdbqt_path
        except Exception as exc:  # noqa: BLE001 - record, never stop the panel
            return {"ligand": ligand_id, "label": entry["label"],
                    "status": f"prep-failed: {exc}"[:200]}
    if pdbqt is None or not Path(pdbqt).is_file():
        return {"ligand": ligand_id, "label": entry["label"],
                "status": "prep-failed: no PDBQT produced"}
    engine = VinaEngine(
        VinaConfig(center=box[0], size=box[1], exhaustiveness=exhaustiveness,
                   num_modes=1, seed=seed, cpu=0, timeout=600.0),
        backend="python", workdir=out_dir / "docking")
    try:
        result = engine.dock(receptor, pdbqt,
                             out_path=out_dir / "docking" / f"{ligand_id}_out.pdbqt")
        if result.best_affinity is not None:
            affinity = result.best_affinity
        else:
            status = f"dock-failed: {result.error}"
    except Exception as exc:  # noqa: BLE001
        status = f"dock-failed: {exc}"[:200]
    return {"ligand": ligand_id, "label": entry["label"], "status": status,
            "best_affinity_kcal_mol": affinity,
            "runtime_s": round(time.perf_counter() - started, 1)}


def read_box(reference_run: Path) -> tuple[tuple, tuple]:
    """Center/size from a previous run's gridbox.txt."""
    values = {}
    for line in (reference_run / "gridbox.txt").read_text().splitlines():
        if "=" in line and not line.startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = float(value)
    center = (values["center_x"], values["center_y"], values["center_z"])
    size = (values["size_x"], values["size_y"], values["size_z"])
    return center, size


def write_report(csv_path: Path, out_dir: Path, panel: list[dict],
                 exhaustiveness: int, seed: int, ratios: tuple[int, int]) -> None:
    """Enrichment metrics + ROC plot from the scored rows."""
    from dockflow_core.enrichment import enrichment_table, plot_roc

    rows = list(load_existing(csv_path).values())
    scored = [row for row in rows
              if row.get("best_affinity_kcal_mol") not in ("", None)]
    if not scored:
        return
    # best affinity per ligand; negated so higher = better ranked
    best: dict[str, float] = {}
    labels: dict[str, bool] = {}
    for row in scored:
        ligand = row["ligand"]
        affinity = float(row["best_affinity_kcal_mol"])
        if ligand not in best or affinity < best[ligand]:
            best[ligand] = affinity
            labels[ligand] = str(row["label"]).lower() in ("true", "1", "yes")
    ligands = sorted(best)
    scores = [-best[name] for name in ligands]
    label_list = [labels[name] for name in ligands]
    metrics = enrichment_table(scores, label_list)
    (out_dir / "enrichment_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    roc_path = plot_roc(scores, label_list, out_dir / "roc.png")

    n_act, n_dec = ratios
    lines = [
        "# DUD-E enrichment benchmark (target hivpr vs 1HVR)",
        "",
        "Protocol: DUD-E hivpr actives + decoys (seeded subsample "
        f"{metrics['n_actives']}/{n_act} actives, "
        f"{metrics['n_decoys']}/{n_dec} decoys of the full 536/35,750 set), "
        "docked with Vina (python bindings) against the 1HVR receptor, "
        "grid box from the XK2 co-crystal ligand (strong assumption), "
        f"exhaustiveness {exhaustiveness}, num_modes 1, seed {seed}.",
        "",
        f"- ligands scored: **{len(scored)}**"
        + (f" ({len(rows) - len(scored)} failed - recorded, not dropped)"
           if len(rows) > len(scored) else ""),
        f"- ROC AUC: **{metrics['roc_auc']:.3f}**",
        f"- EF@1%: {metrics['ef_1pct']} | EF@5%: {metrics['ef_5pct']}",
        f"- BEDROC (alpha 20): {metrics['bedroc_alpha_20']}",
        "",
        "Interpretation: AUC > 0.5 means the Vina score ranks actives "
        "above decoys more often than chance; EF@x% is the actives-fold "
        "enrichment in the top x% of the ranking (1.0 = random). "
        "All metrics are heuristic screening statistics, not binding "
        "free energies.",
        "",
        f"- ROC plot: `{roc_path}`" if roc_path else "- ROC plot: (matplotlib unavailable)",
        f"- per-ligand scores: `{csv_path}`",
        f"- metrics JSON: `{out_dir / 'enrichment_metrics.json'}`",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"metrics: {out_dir / 'enrichment_metrics.json'}")
    print(f"summary: {out_dir / 'summary.md'}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default=str(DEFAULT_DATA))
    parser.add_argument("--out", default=str(Path(__file__).parent / "results"))
    parser.add_argument("--reference-run", default=str(DEFAULT_REF),
                        help="run directory with receptor.pdbqt + gridbox.txt")
    parser.add_argument("--actives", type=int, default=60)
    parser.add_argument("--decoys", type=int, default=500)
    parser.add_argument("--exhaustiveness", type=int, default=4)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--max-seconds", type=float, default=None)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "enrichment_results.csv"

    if args.report_only:
        panel = load_panel(Path(args.data), 10**9, 10**9, args.seed)
        ratios = (sum(1 for p in panel if p["label"]),
                  sum(1 for p in panel if not p["label"]))
        write_report(csv_path, out_dir, panel, args.exhaustiveness,
                     args.seed, ratios)
        return 0

    panel = load_panel(Path(args.data), args.actives, args.decoys, args.seed)
    full = load_panel(Path(args.data), 10**9, 10**9, args.seed)
    ratios = (sum(1 for p in full if p["label"]),
              sum(1 for p in full if not p["label"]))
    existing = load_existing(csv_path)
    todo = [entry for entry in panel if entry["id"] not in existing]
    print(f"enrichment benchmark: {len(panel)} ligands in panel "
          f"({sum(1 for e in panel if e['label'])} actives), "
          f"{len(existing)} scored, {len(todo)} to run | "
          f"exhaustiveness={args.exhaustiveness}, seed={args.seed}", flush=True)

    reference_run = Path(args.reference_run)
    receptor = reference_run / "prepared" / "receptor.pdbqt"
    box = read_box(reference_run)
    started_driver = time.perf_counter()
    for index, entry in enumerate(todo, start=1):
        if args.max_seconds is not None and \
                time.perf_counter() - started_driver > args.max_seconds:
            print(f"time budget reached; {len(todo) - index + 1} ligand(s) "
                  "remain - re-run to resume", flush=True)
            break
        print(f"[{index}/{len(todo)}] {entry['id']} "
              f"({'active' if entry['label'] else 'decoy'}) ...", flush=True)
        row = run_one(entry, out_dir, receptor, box,
                      args.exhaustiveness, args.seed)
        append_row(csv_path, row)
        if row.get("best_affinity_kcal_mol") is not None:
            print(f"    {row['best_affinity_kcal_mol']} kcal/mol "
                  f"({row['runtime_s']}s)", flush=True)
        else:
            print(f"    {row['status']}", flush=True)

    write_report(csv_path, out_dir, panel, args.exhaustiveness,
                 args.seed, ratios)
    print(f"results: {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
