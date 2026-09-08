"""Paper figures from the COMMITTED benchmark artifacts (final pass, item 19).

Nothing here re-runs any experiment: every figure is rendered from the
published results under ``benchmarks/*/results/`` so the paper and the
repository cannot drift apart.  The exact command for each figure is
recorded in ``paper/figures/README.md``.

Outputs: ``paper/figures/fig{1..5}.{png,pdf}``

    python paper/figures/generate_figures.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RED = ROOT / "benchmarks" / "redocking" / "results"
MGL = ROOT / "benchmarks" / "mgltools_comparison" / "results"
ENR = ROOT / "benchmarks" / "enrichment" / "results"
OUT = Path(__file__).parent

SUCCESS = "#1a7f37"
FAIL = "#b35900"
THRESH = "#c62828"


def load_rows(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle) if row.get("pdb_id")]


def _save(fig, name: str) -> None:
    fig.savefig(OUT / f"{name}.png", dpi=200)
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print(f"wrote {OUT / name}.png / .pdf")


def _float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def fig1_redocking_scatter() -> None:
    """Predicted (best-pose) RMSD per complex, coloured by success."""
    rows = load_rows(RED / "redocking_results.csv")
    stats = json.loads((RED / "summary_stats.json").read_text())
    xs, ys, colors, labels = [], [], [], []
    for index, row in enumerate(rows, start=1):
        rmsd = _float(row.get("best_crystal_rmsd_a"))
        if rmsd is None:
            continue
        xs.append(index)
        ys.append(rmsd)
        colors.append(SUCCESS if row.get("success") == "True" else FAIL)
        labels.append(f"{row['pdb_id']}/{row['ligand']}")
    fig, ax = plt.subplots(figsize=(8.5, 5.0), constrained_layout=True)
    ax.vlines(xs, 0, ys, colors=colors, alpha=0.35, lw=1.0, zorder=2)
    ax.scatter(xs, ys, c=colors, s=48, zorder=3)
    for x, y, label in zip(xs, ys, labels, strict=True):
        if y > 2.0:
            ax.annotate(label, (x, y), fontsize=6.5, xytext=(2, 4),
                        textcoords="offset points")
    ax.axhline(2.0, color=THRESH, ls="--", lw=1.2, zorder=2)
    ax.text(len(xs) + 0.3, 2.05, "2.0 Å community-convention threshold",
            fontsize=8, color=THRESH, ha="right")
    ax.set_xlabel("complex index (24-complex set, complexes.csv order)")
    ax.set_ylabel("best-pose symmetry-aware RMSD to crystal (Å)")
    ax.set_ylim(0, max(ys) * 1.12)
    low, high = stats["success_rate_wilson_95ci_pct"]
    ax.set_title(
        f"Redocking benchmark: {stats['successful_within_2a']}/24 complexes "
        f"within 2.0 Å (best-any-rank; 95% CI {low:.1f}–{high:.1f}%)")
    _save(fig, "fig1")


def fig2_rmsd_distribution() -> None:
    """RMSD histogram with the Wilson CI of the success rate."""
    rows = load_rows(RED / "redocking_results.csv")
    stats = json.loads((RED / "summary_stats.json").read_text())
    values = [v for v in (_float(r.get("best_crystal_rmsd_a"))
                          for r in rows) if v is not None]
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    ax.hist(values, bins=12, color="#4c72b0", edgecolor="white", zorder=2)
    ax.axvline(2.0, color=THRESH, ls="--", lw=1.2, zorder=3,
               label="2.0 Å threshold")
    ax.axvline(stats["median_rmsd_a"], color="#55a868", ls="-.", lw=1.4,
               zorder=3, label=f"median {stats['median_rmsd_a']:.2f} Å")
    ax.axvline(stats["mean_rmsd_a"], color="#2f4b7c", ls=":", lw=1.4,
               zorder=3, label=f"mean {stats['mean_rmsd_a']:.2f} Å")
    ax.set_xlabel("best-pose RMSD to crystal (Å)")
    ax.set_ylabel("complexes")
    low, high = stats["success_rate_wilson_95ci_pct"]
    ax.set_title(f"RMSD distribution (n={len(values)}); success "
                 f"{stats['successful_within_2a']}/24 = "
                 f"{stats['success_rate_pct']:.1f}% "
                 f"(95% Wilson CI {low:.1f}–{high:.1f}%)")
    ax.legend(loc="upper right")
    _save(fig, "fig2")


def fig3_mgltools_scatter() -> None:
    """DockFlow vs MGLTools best-pose RMSD, 21 MGL-completed complexes."""
    rows = load_rows(MGL / "comparison.csv")
    pts = []
    for row in rows:
        mgl = _float(row.get("mgl_rmsd_a"))
        df = _float(row.get("df_rmsd_a"))
        if mgl is None or df is None:
            continue
        pts.append((df, mgl, row["pdb_id"]))
    fig, ax = plt.subplots(figsize=(6.0, 6.0), constrained_layout=True)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.scatter(xs, ys, s=55, color="#4c72b0", zorder=3)
    for x, y, pdb_id in pts:
        ax.annotate(pdb_id, (x, y), fontsize=6.5, xytext=(4, 3),
                    textcoords="offset points")
    lim = max(max(xs), max(ys)) * 1.08
    ax.plot([0, lim], [0, lim], color=THRESH, ls="--", lw=1.2,
            label="y = x (identical recovery)", zorder=2)
    ax.set_xlabel("DockFlow arm best-pose RMSD (Å)")
    ax.set_ylabel("MGLTools 1.5.7 arm best-pose RMSD (Å)")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    stats = json.loads((MGL / "comparison_stats.json").read_text())
    ax.set_title(
        f"Preparation-arm comparison, {len(pts)} MGL-completed complexes "
        f"(same box, seed, exhaustiveness)\n"
        f"DockFlow {stats['dockflow_arm_success_within_2a']}/"
        f"{stats['mgl_arm_completed']} vs MGLTools "
        f"{stats['mgl_arm_success_within_2a']}/{stats['mgl_arm_completed']} "
        "within 2.0 Å (shared denominator; 3 MGL crashes excluded)")
    ax.legend(loc="lower right")
    _save(fig, "fig3")


def fig4_hivpr_roc() -> None:
    """DUD-E hivpr ROC curve, re-plotted from the per-ligand scores."""
    metrics = json.loads((ENR / "hivpr" / "enrichment_metrics.json")
                         .read_text())
    rows = list(csv.DictReader(
        open(ENR / "hivpr" / "enrichment_results.csv", newline="",
             encoding="utf-8")))
    scored = [r for r in rows if r.get("best_affinity_kcal_mol")]
    pairs = sorted(
        ((-float(r["best_affinity_kcal_mol"]),
          str(r["label"]).lower() in ("true", "1", "yes"), r["ligand"])
         for r in scored), reverse=True)
    labels = [p[1] for p in pairs]
    n_act = sum(labels)
    n_dec = len(labels) - n_act
    tpr, fpr = [0.0], [0.0]
    tp = fp = 0
    for label in labels:
        if label:
            tp += 1
        else:
            fp += 1
        tpr.append(tp / n_act if n_act else 0.0)
        fpr.append(fp / n_dec if n_dec else 0.0)
    fig, ax = plt.subplots(figsize=(5.5, 5.5), constrained_layout=True)
    ax.plot(fpr, tpr, color="#4c72b0", lw=1.8, zorder=3,
            label=f"DockFlow+Vina (AUC {metrics['roc_auc']:.3f})")
    ax.plot([0, 1], [0, 1], color="#888888", ls="--", lw=1.0,
            label="random (AUC 0.500)", zorder=2)
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_title(f"DUD-E hivpr enrichment evaluation "
                 f"({metrics['n_total']} ligands, exh. 1)")
    ax.legend(loc="lower right")
    _save(fig, "fig4")


def fig5_multitarget_table() -> None:
    """The 3-target enrichment summary as a figure table (item 12)."""
    rows = []
    for target in ("hivpr", "aa2ar", "parp1"):
        path = ENR / target / "enrichment_metrics.json"
        if path.is_file():
            rows.append((target,
                         json.loads(path.read_text(encoding="utf-8"))))
    if not rows:
        print("fig5: no enrichment metrics", file=sys.stderr)
        return
    fig, ax = plt.subplots(figsize=(8.0, 2.4), constrained_layout=True)
    ax.axis("off")
    cell_text = [
        [target, m.get("reference_pdb", ""), str(m["n_actives"]),
         str(m["n_decoys"]), f"{m['roc_auc']:.3f}", str(m["ef_1pct"]),
         str(m["ef_5pct"]), f"{m['bedroc_alpha_20']:.3f}"]
        for target, m in rows
    ]
    table = ax.table(cellText=cell_text,
                     colLabels=["target", "receptor (PDB)", "n actives",
                                "n decoys", "ROC AUC", "EF@1%", "EF@5%",
                                "BEDROC (α=20)"],
                     cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.6)
    ax.set_title("DUD-E enrichment evaluation, 3 targets "
                 "(identical seeded protocol; random AUC = 0.500)",
                 fontsize=10)
    _save(fig, "fig5")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fig1_redocking_scatter()
    fig2_rmsd_distribution()
    fig3_mgltools_scatter()
    fig4_hivpr_roc()
    fig5_multitarget_table()
    return 0


if __name__ == "__main__":
    sys.exit(main())
