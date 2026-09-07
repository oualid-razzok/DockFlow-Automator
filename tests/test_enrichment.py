"""Enrichment metrics tests (audit item 23) with known analytic values."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from dockflow_core.enrichment import (
    bedroc,
    enrichment_factor,
    enrichment_table,
    load_summary_scores,
    plot_roc,
    roc_auc,
)


def test_roc_auc_perfect_and_inverted():
    scores = [10.0, 9.0, 8.0, 2.0, 1.0]
    labels = [1, 1, 1, 0, 0]
    assert roc_auc(scores, labels) == 1.0
    assert roc_auc(scores, [0, 0, 0, 1, 1]) == 0.0


def test_roc_auc_random_expectation():
    """Mean AUC over seeded random rankings converges to 0.5."""
    import random

    aucs = []
    for seed in range(60):
        rng = random.Random(seed)
        labels = [1] * 10 + [0] * 90
        rng.shuffle(labels)
        aucs.append(roc_auc(list(map(float, range(100))), labels))
    assert abs(sum(aucs) / len(aucs) - 0.5) < 0.05


def test_roc_auc_ties():
    # one active tied with two decoys -> 2/3 of comparisons won
    assert roc_auc([1.0, 1.0, 1.0, 0.0], [1, 0, 0, 0]) == pytest.approx(2 / 3)


def test_enrichment_factor_known_values():
    # 10 actives at the top of 100: EF@5% = 5 / (5 * 0.1) = 10
    scores = [float(i) for i in range(100, 0, -1)]
    labels = [i <= 10 for i in range(1, 101)]
    assert enrichment_factor(scores, labels, 0.05) == pytest.approx(10.0)
    # EF@20%: top 20 contains 10 actives -> 10 / (20 * 0.1) = 5
    assert enrichment_factor(scores, labels, 0.20) == pytest.approx(5.0)
    # random ordering -> EF ~ 1
    import random

    rng = random.Random(7)
    shuffled = labels[:]
    rng.shuffle(shuffled)
    ef = enrichment_factor(scores, shuffled, 0.20)
    assert 0.5 < ef < 2.0


def test_bedroc_perfect_and_inverted():
    scores = [float(100 - i) for i in range(100)]
    labels = [i < 10 for i in range(100)]
    assert bedroc(scores, labels, alpha=20.0) == pytest.approx(1.0)
    assert bedroc(scores, [not lab for lab in labels], alpha=20.0) == pytest.approx(0.0)


try:
    import rdkit  # noqa: F401

    _HAS_RDKIT = True
except ImportError:  # pragma: no cover - environments without the prep stack
    _HAS_RDKIT = False


@pytest.mark.skipif(not _HAS_RDKIT, reason="rdkit not installed")
def test_bedroc_matches_rdkit_implementation():
    """Cross-validation against RDKit's reference implementation."""
    import random

    from rdkit.ML.Scoring.Scoring import CalcBEDROC

    for seed in (1, 2, 42):
        rng = random.Random(seed)
        labels = [1] * 10 + [0] * 90
        rng.shuffle(labels)
        # rdkit: (sample, label) sorted best-first, low score = better
        rdkit_input = list(zip(range(100), labels, strict=True))
        expected = CalcBEDROC(rdkit_input, 1, 20.0)
        # ours: higher score = better
        ours = bedroc([-float(i) for i in range(100)], labels, alpha=20.0)
        assert ours == pytest.approx(expected, abs=1e-12)


def test_bedroc_alpha_validation():
    with pytest.raises(ValueError, match="alpha"):
        bedroc([1.0, 2.0], [1, 0], alpha=0.0)


def test_metrics_reject_degenerate_inputs():
    with pytest.raises(ValueError, match="actives"):
        roc_auc([1.0, 2.0], [0, 0])
    with pytest.raises(ValueError, match="decoys"):
        bedroc([1.0, 2.0], [1, 1])
    with pytest.raises(ValueError, match="same length"):
        enrichment_factor([1.0], [1, 0], 0.1)


def test_enrichment_table_shape():
    table = enrichment_table([10.0, 9.0, 1.0, 0.5], [1, 1, 0, 0])
    assert table["n_actives"] == 2 and table["n_decoys"] == 2
    assert table["roc_auc"] == 1.0
    assert table["ef_1pct"] == 2.0  # top 1% of 4 = 1 molecule, 1 active
    assert "bedroc_alpha_20" in table


def test_load_summary_scores(tmp_path: Path):
    csv_path = tmp_path / "summary.csv"
    csv_path.write_text(
        "ligand,pose,affinity_kcal_mol,rmsd_lb,rmsd_ub,crystal_rmsd,"
        "docking_score_efficiency,num_heavy_atoms,runtime_s,backend,error\n"
        "xk2,1,-9.4,0,0,1.2,-0.2,47,10,python,\n"
        "xk2,2,-8.7,1.2,2.0,, ,47,10,python,\n"
        "aspirin,1,-5.1,0,0,, ,13,5,python,\n"
        "caffeine,1,-4.9,0,0,, ,24,4,python,\n",
        encoding="utf-8",
    )
    scores, labels, ligands = load_summary_scores(csv_path, ["xk2"])
    assert ligands == ["aspirin", "caffeine", "xk2"]
    assert labels == [False, False, True]
    assert scores[2] == 9.4  # best (most negative) affinity, negated
    with pytest.raises(Exception, match="none of the actives"):
        load_summary_scores(csv_path, ["not_there"])


def test_plot_roc_optional(tmp_path: Path):
    path = plot_roc([3.0, 2.0, 1.0, 0.5], [1, 1, 0, 0], tmp_path / "roc.png")
    if path is not None:  # matplotlib installed
        assert path.is_file() and path.stat().st_size > 1000


def test_bedroc_random_baseline_documented():
    """The random BEDROC at alpha=20 is NOT 0.5 (only true as alpha->0).

    Documents the behaviour the docstring warns about; guards against an
    accidental 'normalisation' that would silently change the metric.
    """
    import random

    values = []
    for seed in range(20):
        rng = random.Random(seed)
        labels = [1] * 10 + [0] * 90
        rng.shuffle(labels)
        values.append(bedroc([-float(i) for i in range(100)], labels, alpha=20.0))
    mean = sum(values) / len(values)
    assert not math.isclose(mean, 0.5, abs_tol=0.1)
    assert 0.0 < mean < 0.3
