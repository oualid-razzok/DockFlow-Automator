"""Virtual-screening enrichment metrics (audit item 23).

Implements the standard enrichment metrics used to quantify whether known
actives score better than decoys:

* **ROC AUC** - area under the receiver-operating curve (probability that a
  random active scores better than a random decoy; 0.5 = random).
* **EF@x%** - enrichment factor in the top x% of the ranked list
  (actives fraction in the top x% divided by the overall actives fraction).
* **BEDROC** (alpha=20) - Boltzmann-enhanced discrimination of ROC, which
  weights early recognition exponentially; 1 = perfect, 0 = inverted.
  Implemented exactly as defined by Truchon & Bayly (2007), JCIM 47,
  488-508, and verified to reproduce RDKit's ``ML.Scoring.Scoring``
  implementation to machine precision.  Note the random baseline is 0.5
  only in the alpha->0 limit; at alpha=20 it depends on the actives
  fraction (e.g. ~0.15 at 10% actives), so compare BEDROC values between
  rankings with the same alpha and actives count, or against a
  label-shuffled null:

  ::

      RIE     = sum(exp(-alpha * rank_i / N) for active i) / (n_act * denom)
      denom   = (1 - exp(-alpha)) / (N * (exp(alpha/N) - 1))
      RIEmax  = (1 - exp(-alpha * n_act/N)) / ((n_act/N) * (1 - exp(-alpha)))
      RIEmin  = (1 - exp( alpha * n_act/N)) / ((n_act/N) * (1 - exp( alpha)))
      BEDROC  = (RIE - RIEmin) / (RIEmax - RIEmin)

Scores are always "higher is better" (e.g. negated Vina affinities); ranks
are 1-based with rank 1 = best.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from pathlib import Path

from .utils import DockFlowError

__all__ = [
    "roc_auc",
    "enrichment_factor",
    "bedroc",
    "enrichment_table",
    "load_summary_scores",
    "plot_roc",
]


def _validate(scores: Sequence[float], labels: Sequence[bool | int]) -> None:
    if len(scores) != len(labels):
        raise ValueError("scores and labels must have the same length")
    if not scores:
        raise ValueError("no molecules given")
    n_actives = sum(1 for label in labels if label)
    if n_actives == 0:
        raise ValueError("no actives given; enrichment is undefined")
    if n_actives == len(labels):
        raise ValueError("no decoys given; enrichment is undefined")


def _ranked(scores: Sequence[float], labels: Sequence[bool | int]) -> list[bool]:
    """Labels sorted by score descending (best first), stable on ties."""
    pairs = sorted(
        zip(scores, labels, strict=True),
        key=lambda pair: -float(pair[0]),
    )
    return [bool(label) for _, label in pairs]


def roc_auc(scores: Sequence[float], labels: Sequence[bool | int]) -> float:
    """Area under the ROC curve via the rank-sum (Mann-Whitney) statistic."""
    _validate(scores, labels)
    n = len(scores)
    n_actives = sum(1 for label in labels if label)
    n_decoys = n - n_actives
    # Rank ascending (average ranks for ties) - classic AUC formulation.
    order = sorted(range(n), key=lambda i: float(scores[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and float(scores[order[j + 1]]) == float(scores[order[i]]):
            j += 1
        average_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = average_rank
        i = j + 1
    rank_sum = sum(rank for rank, label in zip(ranks, labels, strict=True) if label)
    return float((rank_sum - n_actives * (n_actives + 1) / 2.0) / (n_actives * n_decoys))


def enrichment_factor(
    scores: Sequence[float], labels: Sequence[bool | int], fraction: float
) -> float:
    """Enrichment factor in the top ``fraction`` of the ranked list.

    EF@x% = (actives within top x%) / (x% * N) divided by (n_actives / N).
    Returns the maximum representable value (n_act / min(n_act, top_n))
    when the top slice contains every active.
    """
    _validate(scores, labels)
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"fraction must be in (0, 1], got {fraction}")
    ranked_labels = _ranked(scores, labels)
    n = len(ranked_labels)
    n_actives = sum(ranked_labels)
    top_n = max(1, math.ceil(n * fraction))
    actives_in_top = sum(ranked_labels[:top_n])
    expected_in_top = top_n * (n_actives / n)
    return float(actives_in_top / expected_in_top)


def bedroc(
    scores: Sequence[float], labels: Sequence[bool | int], alpha: float = 20.0
) -> float:
    """Boltzmann-Enhanced Discrimination of ROC (Truchon & Bayly 2007)."""
    _validate(scores, labels)
    if alpha <= 0.0:
        raise ValueError("alpha must be > 0")
    ranked_labels = _ranked(scores, labels)
    n = len(ranked_labels)
    n_actives = sum(ranked_labels)
    ratio = n_actives / n
    denom = (1.0 - math.exp(-alpha)) / (n * (math.exp(alpha / n) - 1.0))
    sum_exp = sum(
        math.exp(-(alpha * (rank + 1)) / n)
        for rank, label in enumerate(ranked_labels)
        if label
    )
    rie = sum_exp / (n_actives * denom)
    rie_max = (1.0 - math.exp(-alpha * ratio)) / (ratio * (1.0 - math.exp(-alpha)))
    rie_min = (1.0 - math.exp(alpha * ratio)) / (ratio * (1.0 - math.exp(alpha)))
    if rie_max == rie_min:
        return 1.0
    return float((rie - rie_min) / (rie_max - rie_min))


def enrichment_table(
    scores: Sequence[float], labels: Sequence[bool | int], alpha: float = 20.0
) -> dict[str, float | int]:
    """All enrichment metrics for one actives+decoys ranking."""
    _validate(scores, labels)
    return {
        "n_total": len(scores),
        "n_actives": sum(1 for label in labels if label),
        "n_decoys": sum(1 for label in labels if not label),
        "roc_auc": round(roc_auc(scores, labels), 4),
        "ef_1pct": round(enrichment_factor(scores, labels, 0.01), 3),
        "ef_5pct": round(enrichment_factor(scores, labels, 0.05), 3),
        f"bedroc_alpha_{alpha:g}": round(bedroc(scores, labels, alpha=alpha), 4),
    }


# ---------------------------------------------------------------------------
# summary.csv loading
# ---------------------------------------------------------------------------
def load_summary_scores(
    summary_csv: str | Path,
    actives: Sequence[str] | set[str],
) -> tuple[list[float], list[bool], list[str]]:
    """Best (most negative) affinity per ligand from a DockFlow summary.csv.

    Returns ``(scores, labels, ligands)`` where scores are **negated**
    affinities so that "higher is better" and labels mark actives.
    """
    actives_set = {str(active).strip().lower() for active in actives}
    best: dict[str, float] = {}
    with open(summary_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"ligand", "affinity_kcal_mol"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise DockFlowError(
                f"{summary_csv} is not a DockFlow summary.csv "
                f"(columns: {reader.fieldnames})"
            )
        for row in reader:
            ligand = (row.get("ligand") or "").strip()
            raw = (row.get("affinity_kcal_mol") or "").strip()
            if not ligand or not raw:
                continue
            try:
                affinity = float(raw)
            except ValueError:
                continue
            if ligand not in best or affinity < best[ligand]:
                best[ligand] = affinity
    if not best:
        raise DockFlowError(f"no scored ligands found in {summary_csv}")
    unknown = [name for name in actives_set if name not in
               {key.lower() for key in best}]
    ligands = sorted(best)
    scores = [-best[name] for name in ligands]
    labels = [name.lower() in actives_set for name in ligands]
    if not any(labels):
        raise DockFlowError(
            "none of the actives appear in the summary (checked "
            f"{len(actives_set)} ids; unknown: {', '.join(sorted(unknown)[:10])})"
        )
    return scores, labels, ligands


def plot_roc(
    scores: Sequence[float], labels: Sequence[bool | int], out_path: str | Path
) -> Path | None:
    """Render the ROC curve to a PNG (matplotlib optional)."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    _validate(scores, labels)
    pairs = sorted(zip(scores, labels, strict=True), key=lambda pair: -float(pair[0]))
    n_actives = sum(1 for _, label in pairs if label)
    n_decoys = len(pairs) - n_actives
    x, y = [0.0], [0.0]
    tp = fp = 0
    for _, label in pairs:
        if label:
            tp += 1
        else:
            fp += 1
        x.append(fp / n_decoys)
        y.append(tp / n_actives)
    fig, axis = plt.subplots(figsize=(5.2, 5.2), constrained_layout=True)
    axis.plot(x, y, color="#1a7f37", lw=2, label="ranking")
    axis.plot([0, 1], [0, 1], "--", color="#888888", lw=1, label="random")
    axis.set_xlabel("false positive rate (decoys)")
    axis.set_ylabel("true positive rate (actives)")
    axis.set_title(f"ROC (AUC = {roc_auc(scores, labels):.3f})")
    axis.legend(loc="lower right")
    axis.set_xlim(-0.02, 1.02)
    axis.set_ylim(-0.02, 1.02)
    out = Path(out_path)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out
