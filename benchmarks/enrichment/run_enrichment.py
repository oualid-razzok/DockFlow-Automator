"""DUD-E enrichment benchmark runner (work order item 3; final pass items 10-12).

Docks a DUD-E actives + decoys panel against each target's reference
receptor with a ligand-derived grid box and publishes the enrichment
metrics (ROC AUC, EF@1%, EF@5%, BEDROC).

Targets (final pass, item 12 - at least 3 for a small panel):

* ``hivpr`` - receptor 1HVR (the DockFlow example), box from XK2;
* ``aa2ar`` - receptor 3EML, box from the ZMA co-crystal ligand;
* ``parp1`` - receptor 2RD6, box from the 78P co-crystal ligand.

Protocol (identical for every target, documented for reproduction in
``README.md``):

* panel: seeded subsample (``random.Random(seed).sample``) of 40 actives
  and 93 decoys from the DUD-E ``actives_final.ism`` /
  ``decoys_final.ism`` lists (the ISM files are committed under ``data/``
  with their SHA-256 recorded in README.md, so the exact panel is
  reproducible);
* ligand preparation: Meeko (embed 3D + MMFF minimisation, seed),
  input protonation state, no tautomer enumeration, salt stripping by
  Meeko default;
* receptor preparation: DockFlow (engine auto -> openbabel here),
  Gasteiger charges, deposited protonation + template hydrogens;
* grid box: from the reference structure's co-crystallised ligand,
  4.0 A padding (strong assumption);
* docking: Vina python bindings, exhaustiveness 1, num_modes 1,
  seed 2026 (screening-grade settings);
* scoring: the single Vina ``affinity`` column (lower = better);
  score ties are broken by ascending ligand identifier (deterministic).

The runner is resumable like the redocking benchmark:

* one CSV row per ligand, appended the moment it is scored,
* prepared ligand PDBQTs are cached (re-runs skip them),
* ``--report-only`` recomputes metrics/plots from the CSV,
* ``--max-seconds`` stops cleanly between ligands for chunked execution,
* ``--reproduce`` (hivpr only) re-docks the exact 133-ligand panel in a
  scratch directory and asserts the metrics match the published ones
  within 1e-6.

Usage::

    python benchmarks/enrichment/run_enrichment.py --target hivpr
    python benchmarks/enrichment/run_enrichment.py --target aa2ar
    python benchmarks/enrichment/run_enrichment.py --target parp1
    python benchmarks/enrichment/run_enrichment.py --reproduce
    python benchmarks/enrichment/run_enrichment.py --report-only --target hivpr
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

DEFAULT_DATA = Path(__file__).parent / "data"
DEFAULT_OUT_ROOT = Path(__file__).parent / "results"
DEFAULT_REF_RUN = (Path(__file__).resolve().parents[1] / "redocking"
                   / "results" / "1HVR_benchmark")

DUDE_BASE = "https://dude.docking.org/targets/{target}"

# Per-target protocol parameters (final pass, item 12).  The reference
# PDB + co-crystal ligand define the receptor and the grid box exactly
# like the hivpr protocol (box from the co-crystal ligand = strong
# assumption).  DUD-E full-set sizes are recorded for the subsample
# ratio in the summary.
TARGETS = {
    "hivpr": {
        "reference_pdb": "1HVR",
        "reference_ligand": "XK2",
        "dude_full_actives": 536,
        "dude_full_decoys": 35750,
        "note": "DockFlow example target; receptor + box reused from the "
                "redocking benchmark run (identical preparation).",
    },
    "aa2ar": {
        "reference_pdb": "3EML",
        "reference_ligand": "ZMA",
        "dude_full_actives": 482,
        "dude_full_decoys": 31550,
        "note": "adenosine A2A receptor; reference structure 3EML "
                "(ZM241385 co-crystal).",
    },
    "parp1": {
        "reference_pdb": "2RD6",
        "reference_ligand": "78P",
        "dude_full_actives": 508,
        "dude_full_decoys": 30050,
        "note": "poly(ADP-ribose) polymerase 1; reference structure 2RD6 "
                "(inhibitor 78P co-crystal).",
    },
}

FIELDNAMES = ["ligand", "label", "status", "best_affinity_kcal_mol",
              "runtime_s"]

N_ACTIVES = 40   # seeded subsample size (matches the published hivpr run)
N_DECOYS = 93    # ditto


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def ensure_dude_data(data_dir: Path, target: str) -> tuple[Path, Path]:
    """Committed (or downloaded) DUD-E ISM lists for one target."""
    actives = data_dir / f"{target}_actives_final.ism"
    decoys = data_dir / f"{target}_decoys_final.ism"
    import requests

    for path, kind in ((actives, "actives"), (decoys, "decoys")):
        if path.is_file():
            continue
        url = f"{DUDE_BASE.format(target=target)}/{kind}_final.ism"
        print(f"downloading {url} ...", flush=True)
        response = requests.get(url, timeout=180)
        response.raise_for_status()
        data_dir.mkdir(parents=True, exist_ok=True)
        path.write_text(response.text, encoding="utf-8")
    return actives, decoys


def load_panel(data_dir: Path, target: str, n_actives: int, n_decoys: int,
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

    actives = parse(data_dir / f"{target}_actives_final.ism", True)
    decoys = parse(data_dir / f"{target}_decoys_final.ism", False)
    if n_actives < len(actives):
        actives = rng.sample(actives, n_actives)
    if n_decoys < len(decoys):
        decoys = rng.sample(decoys, n_decoys)
    return actives + decoys


def duplicate_report(panel: list[dict]) -> dict:
    """Duplicate-SMILES report (item 10: duplicate handling documented).

    DUD-E distributes already-deduplicated lists; this check makes that
    claim verifiable instead of trusted.  Duplicates are REPORTED, never
    dropped - dropping would silently change the published panel.
    Canonicalisation uses RDKit when importable, exact SMILES otherwise.
    """
    try:
        from rdkit import Chem

        keys = [Chem.CanonSmiles(entry["smiles"]) for entry in panel]
        method = "rdkit-canonical-smiles"
    except ImportError:
        keys = [entry["smiles"] for entry in panel]
        method = "exact-smiles-string"
    counts: dict[str, int] = {}
    for key in keys:
        counts[key] = counts.get(key, 0) + 1
    duplicates = {key: n for key, n in counts.items() if n > 1}
    return {
        "dedup_method": method,
        "duplicate_smiles": duplicates,
        "n_duplicate_entries": sum(n - 1 for n in duplicates.values()),
    }


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


# ---------------------------------------------------------------------------
# receptor + grid box per target
# ---------------------------------------------------------------------------
def prepare_receptor_and_box(target: str, out_dir: Path,
                             reference_run: Path | None) -> tuple[Path, tuple]:
    """(receptor PDBQT, (center, size)) for one target.

    ``reference_run`` (the hivpr protocol): reuse the redocking
    benchmark's prepared receptor + grid box verbatim - the enrichment
    run must be bit-identical to the published one.  Otherwise (aa2ar,
    parp1): download the reference PDB, prepare the receptor exactly
    like the redocking benchmark (engine auto, Gasteiger), derive the
    box from the co-crystallised ligand with 4.0 A padding.
    """
    prepared = out_dir / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    receptor = prepared / "receptor.pdbqt"
    config_path = out_dir / "gridbox.txt"

    if reference_run is not None and \
            (reference_run / "prepared" / "receptor.pdbqt").is_file() \
            and (reference_run / "gridbox.txt").is_file():
        if not receptor.is_file():
            receptor.write_bytes(
                (reference_run / "prepared" / "receptor.pdbqt").read_bytes())
        if not config_path.is_file():
            config_path.write_text(
                (reference_run / "gridbox.txt").read_text(encoding="utf-8"),
                encoding="utf-8")
        return receptor, read_box(reference_run)

    spec = TARGETS[target]
    pdb_id = spec["reference_pdb"]
    resname = spec["reference_ligand"]
    raw = out_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    pdb_path = raw / f"{pdb_id.lower()}.pdb"
    if not pdb_path.is_file():
        import requests

        response = requests.get(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=120)
        response.raise_for_status()
        pdb_path.write_text(response.text, encoding="utf-8")

    if not receptor.is_file():
        from dockflow_core.preparator import (
            ReceptorPreparator,
            ReceptorPrepOptions,
        )

        result = ReceptorPreparator(ReceptorPrepOptions(
            engine="auto", charge_model="gasteiger")).prepare(
                pdb_path, prepared, basename="receptor")
        receptor = result.pdbqt_path

    # grid box: co-crystallised ligand + 4.0 A padding (strong assumption)
    if not config_path.is_file():
        from dockflow_core.analyzer import extract_reference_atoms
        from dockflow_core.gridbox import GridBox, box_from_atoms

        ligand_atoms = extract_reference_atoms(pdb_path, resname)
        box: GridBox = box_from_atoms(ligand_atoms, padding=4.0,
                                      source=f"ligand:{resname}")
        lines = [
            f"# DUD-E enrichment grid box - target {target} - "
            f"reference {pdb_id} ligand {resname} (4.0 A padding, "
            "assumption: strong)",
            f"center_x = {box.center[0]}",
            f"center_y = {box.center[1]}",
            f"center_z = {box.center[2]}",
            f"size_x = {box.size[0]}",
            f"size_y = {box.size[1]}",
            f"size_z = {box.size[2]}",
            "exhaustiveness = 1",
        ]
        config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return receptor, read_box(out_dir)


def read_box(reference_dir: Path) -> tuple[tuple, tuple]:
    """Center/size from a run's gridbox.txt."""
    values = {}
    for line in (reference_dir / "gridbox.txt").read_text().splitlines():
        if "=" in line and not line.startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = float(value)
    center = (values["center_x"], values["center_y"], values["center_z"])
    size = (values["size_x"], values["size_y"], values["size_z"])
    return center, size


# ---------------------------------------------------------------------------
# docking + reporting
# ---------------------------------------------------------------------------
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


def _random_baselines(scores: list[float], labels: list[bool],
                      alpha: float = 20.0, n_permutations: int = 10000
                      ) -> dict:
    """Empirical random baselines (item 10: printed explicitly).

    The BEDROC random baseline is NOT 0.5 at alpha = 20 - it depends on
    the actives fraction.  Instead of trusting a formula, the baseline
    is measured: mean metric over ``n_permutations`` seeded random
    shuffles of the same labels.  Deterministic (seed fixed).
    """
    from dockflow_core.enrichment import bedroc, roc_auc

    rng = random.Random(2026)
    auc_values = []
    bedroc_values = []
    labels_array = list(labels)
    for _ in range(n_permutations):
        rng.shuffle(labels_array)
        auc_values.append(roc_auc(scores, labels_array))
        bedroc_values.append(bedroc(scores, labels_array, alpha=alpha))
    return {
        "roc_auc_random": round(sum(auc_values) / len(auc_values), 4),
        f"bedroc_alpha_{alpha:g}_random": round(
            sum(bedroc_values) / len(bedroc_values), 4),
        "ef_random": 1.0,
        "n_permutations": n_permutations,
    }


def compute_metrics(csv_path: Path) -> tuple[dict, list[str], list[bool]] | None:
    """(metrics, ligand ids, labels) from the scored CSV rows."""
    from dockflow_core.enrichment import enrichment_table

    rows = list(load_existing(csv_path).values())
    scored = [row for row in rows
              if row.get("best_affinity_kcal_mol") not in ("", None)]
    if not scored:
        return None
    best: dict[str, float] = {}
    labels: dict[str, bool] = {}
    for row in scored:
        ligand = row["ligand"]
        affinity = float(row["best_affinity_kcal_mol"])
        if ligand not in best or affinity < best[ligand]:
            best[ligand] = affinity
            labels[ligand] = str(row["label"]).lower() in ("true", "1", "yes")
    ligands = sorted(best)  # ties in score break by ascending ligand id
    scores = [-best[name] for name in ligands]  # negated: higher = better
    label_list = [labels[name] for name in ligands]
    return enrichment_table(scores, label_list), ligands, label_list


def write_report(csv_path: Path, out_dir: Path, target: str, panel: list[dict],
                 exhaustiveness: int, seed: int) -> dict | None:
    """Enrichment metrics + ROC plot from the scored rows (items 10, 11)."""
    from dockflow_core.enrichment import plot_roc

    result = compute_metrics(csv_path)
    if result is None:
        return None
    metrics, ligands, label_list = result
    spec = TARGETS[target]
    n_act, n_dec = spec["dude_full_actives"], spec["dude_full_decoys"]
    dupes = duplicate_report(panel)
    # negated affinities so higher = better ranked; ties break by
    # ascending ligand id (sorted) - deterministic, documented
    scores = [-float(load_existing(csv_path)[name]["best_affinity_kcal_mol"])
              for name in ligands]
    baselines = _random_baselines(scores, label_list)
    metrics.update({
        "target": target,
        "reference_pdb": spec["reference_pdb"],
        "reference_ligand": spec["reference_ligand"],
        "exhaustiveness": exhaustiveness,
        "seed": seed,
        "subsample": {
            "actives": f"{metrics['n_actives']}/{n_act}",
            "decoys": f"{metrics['n_decoys']}/{n_dec}",
            "rng": "random.Random(seed).sample",
            "seed": seed,
        },
        "duplicates": dupes,
        "random_baselines": baselines,
    })
    (out_dir / "enrichment_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    roc_path = plot_roc(scores, label_list, out_dir / "roc.png")

    n_total = metrics["n_actives"] + metrics["n_decoys"]
    lines = [
        f"# DUD-E enrichment evaluation (target {target}, subsample "
        f"{n_total} ligands)",
        "",
        f"Enrichment evaluation on the DUD-E {target} subsample "
        f"({n_total} ligands): ROC AUC {metrics['roc_auc']:.3f}.  This is "
        "an **EVALUATION** of the DockFlow+Vina pipeline on one target, "
        "NOT a claim of superior virtual-screening performance "
        "(single-target, exhaustiveness 1, no comparative baseline against "
        "other tools on the same panel).",
        "",
        f"Protocol: DUD-E {target} actives + decoys (seeded subsample "
        f"{metrics['n_actives']}/{n_act} actives, {metrics['n_decoys']}/{n_dec} "
        f"decoys of the full set), receptor {spec['reference_pdb']} prepared "
        "by DockFlow (engine auto), grid box from the "
        f"{spec['reference_ligand']} co-crystal ligand (strong assumption, "
        f"4.0 A padding), Vina python bindings, exhaustiveness "
        f"{exhaustiveness}, num_modes 1, seed {seed}.  Ligand preparation: "
        "Meeko (embed 3D + MMFF minimise, seeded), input protonation "
        "states, no tautomer enumeration, Meeko-default salt stripping.  "
        f"{spec['note']}",
        "",
        f"- ligands scored: **{n_total}**",
        f"- ROC AUC: **{metrics['roc_auc']:.3f}** "
        f"(random baseline {baselines['roc_auc_random']:.3f})",
        f"- EF@1%: {metrics['ef_1pct']} | EF@5%: {metrics['ef_5pct']} "
        "(random baseline 1.0)",
        f"- BEDROC (alpha 20): {metrics['bedroc_alpha_20']} "
        f"(random baseline {baselines['bedroc_alpha_20_random']} - the "
        "BEDROC random baseline is NOT 0.5 at alpha 20; it depends on the "
        "actives fraction, so it is measured here over "
        f"{baselines['n_permutations']} seeded random permutations of the "
        "same panel)",
        "",
        "## Duplicate handling (item 10)",
        "",
        f"- DUD-E distributes already-deduplicated lists; this runner "
        f"VERIFIES that (method: {dupes['dedup_method']}) instead of "
        "trusting it.  Duplicates are reported, never dropped - the "
        "published panel is exactly the seeded subsample.",
        f"- duplicate entries found: {dupes['n_duplicate_entries']}",
        "",
        "## Interpretation",
        "",
        "AUC > 0.5 means the Vina score ranks actives above decoys more "
        "often than chance; EF@x% is the actives-fold enrichment in the "
        "top x% of the ranking (1.0 = random).  All metrics are heuristic "
        "screening statistics, not binding free energies.  Caveats: "
        "single exhaustiveness-1 run, single target, no comparative "
        "baseline - the AUC is informative for 'does the pipeline enrich "
        "actives at all', not for 'is DockFlow better than tool X'.  "
        "Three targets is still a small panel; a publication-grade VS "
        "benchmark would use the full DUD-E (102 targets) - see "
        "docs/SCIENTIFIC_VALIDATION.md.",
        "",
        f"- ROC plot: `{roc_path}`" if roc_path
        else "- ROC plot: (matplotlib unavailable)",
        f"- per-ligand scores: `{csv_path}`",
        f"- metrics JSON: `{out_dir / 'enrichment_metrics.json'}`",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    print(f"metrics: {out_dir / 'enrichment_metrics.json'}")
    print(f"summary: {out_dir / 'summary.md'}")
    return metrics


def write_multi_target_summary(out_root: Path) -> None:
    """The 3-target comparison table (item 12)."""
    rows = []
    for target in TARGETS:
        metrics_path = out_root / target / "enrichment_metrics.json"
        if not metrics_path.is_file():
            continue
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        rows.append((target, metrics))
    if not rows:
        return
    lines = [
        "# DUD-E enrichment evaluation - multi-target summary (item 12)",
        "",
        "| target | receptor (PDB) | n actives | n decoys | ROC AUC | "
        "EF@1% | EF@5% | BEDROC (a=20) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for target, metrics in rows:
        lines.append(
            f"| {target} | {metrics.get('reference_pdb')} "
            f"| {metrics['n_actives']} | {metrics['n_decoys']} "
            f"| {metrics['roc_auc']:.3f} | {metrics['ef_1pct']} "
            f"| {metrics['ef_5pct']} | {metrics['bedroc_alpha_20']} |")
    lines += [
        "",
        "All targets share the identical protocol (see each target's "
        "summary.md): seeded 40/93 subsample, Meeko ligand prep, "
        "ligand-derived grid box (4.0 A padding), Vina exhaustiveness 1, "
        "num_modes 1, seed 2026.  ROC AUC random baseline 0.500; "
        "EF random baseline 1.0; the BEDROC random baseline depends on "
        "the actives fraction (recorded per target in "
        "enrichment_metrics.json).",
        "",
        "**Limitation (explicit):** three targets is still a small panel "
        "and every number above is a single exhaustiveness-1 seeded run.  "
        "A publication-grade virtual-screening benchmark would evaluate "
        "the full DUD-E (102 targets) with multiple seeds - that is "
        "future work (see ROADMAP_POST_V1.md), stated here so the table "
        "is not over-read.",
        "",
    ]
    (out_root / "summary.md").write_text("\n".join(lines) + "\n",
                                         encoding="utf-8")
    print(f"multi-target summary: {out_root / 'summary.md'}")


# ---------------------------------------------------------------------------
# reproduce mode (item 10)
# ---------------------------------------------------------------------------
def reproduce(data_dir: Path, out_root: Path, reference_run: Path,
              exhaustiveness: int, seed: int) -> int:
    """Re-run the exact hivpr panel and assert the published metrics.

    Docks the exact 133-ligand panel into a scratch directory
    (``results-reproduce/``, resumable like every other runner: ligands
    already scored in its CSV are skipped and prepared PDBQTs are
    cached) and compares against the committed
    ``enrichment_metrics.json``.  Vina's
    search is thread-count-dependent (same lesson as the v0.3.1 CI fix):
    on the same machine class (2 vCPU, cpu=0) the metrics reproduce to
    1e-6; on a different thread count the search trajectories differ and
    the assertion will say so honestly instead of silently passing.
    """
    scratch = out_root.parent / "results-reproduce"
    scratch.mkdir(parents=True, exist_ok=True)
    target = "hivpr"
    panel = load_panel(data_dir, target, N_ACTIVES, N_DECOYS, seed)
    receptor, box = prepare_receptor_and_box(target, scratch, reference_run)
    csv_path = scratch / "enrichment_results.csv"
    existing = load_existing(csv_path)
    todo = [entry for entry in panel
            if entry["id"] not in existing]
    print(f"reproduce: docking the exact {len(panel)}-ligand panel into "
          f"{scratch} ({len(existing)} already scored, {len(todo)} to run; "
          f"exhaustiveness {exhaustiveness}, seed {seed}) ...", flush=True)
    for index, entry in enumerate(todo, start=1):
        print(f"[{index}/{len(todo)}] {entry['id']} "
              f"({'active' if entry['label'] else 'decoy'})", flush=True)
        append_row(csv_path, run_one(entry, scratch, receptor, box,
                                     exhaustiveness, seed))
    fresh = compute_metrics(csv_path)[0]
    published_path = out_root / target / "enrichment_metrics.json"
    if not published_path.is_file():
        print(f"FAIL: no published metrics at {published_path}",
              file=sys.stderr)
        return 1
    published = json.loads(published_path.read_text(encoding="utf-8"))
    keys = ["roc_auc", "ef_1pct", "ef_5pct", "bedroc_alpha_20"]
    mismatches = [
        (key, published.get(key), fresh.get(key))
        for key in keys
        if abs(float(published.get(key, 0.0)) - float(fresh.get(key, 0.0)))
        > 1e-6
    ]
    if mismatches:
        print("REPRODUCE FAILED - metrics differ beyond 1e-6:", file=sys.stderr)
        for key, want, got in mismatches:
            print(f"  {key}: published {want} vs fresh {got}",
                  file=sys.stderr)
        return 1
    print("REPRODUCE PASSED: all metrics match the published values "
          "within 1e-6")
    return 0


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", default="hivpr",
                        choices=sorted(TARGETS),
                        help="DUD-E target to evaluate")
    parser.add_argument("--data", default=str(DEFAULT_DATA))
    parser.add_argument("--out-root", default=str(DEFAULT_OUT_ROOT))
    parser.add_argument("--reference-run", default=str(DEFAULT_REF_RUN),
                        help="redocking run dir whose receptor + grid box "
                             "to reuse (hivpr protocol)")
    parser.add_argument("--actives", type=int, default=N_ACTIVES)
    parser.add_argument("--decoys", type=int, default=N_DECOYS)
    parser.add_argument("--exhaustiveness", type=int, default=1)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--max-seconds", type=float, default=None)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--reproduce", action="store_true",
                        help="re-run the exact published hivpr panel and "
                             "assert the metrics match within 1e-6")
    args = parser.parse_args()

    data_dir = Path(args.data)
    out_root = Path(args.out_root)
    target = args.target

    if args.reproduce:
        return reproduce(data_dir, out_root, Path(args.reference_run),
                         args.exhaustiveness, args.seed)

    out_dir = out_root / target
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "enrichment_results.csv"

    if args.report_only:
        panel = load_panel(data_dir, target, 10**9, 10**9, args.seed)
        if write_report(csv_path, out_dir, target, panel,
                        args.exhaustiveness, args.seed) is None:
            print("no scored rows yet", file=sys.stderr)
            return 1
        write_multi_target_summary(out_root)
        return 0

    ensure_dude_data(data_dir, target)
    panel = load_panel(data_dir, target, args.actives, args.decoys, args.seed)
    reference_run = (Path(args.reference_run) if target == "hivpr"
                     else None)
    receptor, box = prepare_receptor_and_box(target, out_dir, reference_run)
    existing = load_existing(csv_path)
    todo = [entry for entry in panel if entry["id"] not in existing]
    print(f"enrichment benchmark ({target}): {len(panel)} ligands in panel "
          f"({sum(1 for e in panel if e['label'])} actives), "
          f"{len(existing)} scored, {len(todo)} to run | "
          f"exhaustiveness={args.exhaustiveness}, seed={args.seed}", flush=True)

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

    write_report(csv_path, out_dir, target, panel,
                 args.exhaustiveness, args.seed)
    write_multi_target_summary(out_root)
    print(f"results: {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
