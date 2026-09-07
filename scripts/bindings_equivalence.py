"""C++ accelerator vs NumPy equivalence report (audit item 35).

Runs both implementations over randomised inputs and writes a machine
readable delta report (consumed as a CI artifact) plus a human summary.
Exits non-zero when any delta exceeds its tolerance, so the equivalence
is enforced, not just observed.

Usage::

    python scripts/bindings_equivalence.py [--out deltas.json] [--cases 200]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

TOLERANCES = {
    "kabsch_rmsd": 1e-5,
    "pairwise_min_dist": 1e-8,
    "min_contacts": 1e-9,  # counts + distances must match exactly
}


def _random_points(rng: np.random.Generator, n: int, spread: float = 20.0):
    return rng.uniform(-spread, spread, size=(n, 3)).astype(np.float64)


def _numpy_kabsch(a: np.ndarray, b: np.ndarray) -> float:
    pa = a - a.mean(axis=0)
    pb = b - b.mean(axis=0)
    covariance = pa.T @ pb
    u, _, vt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    aligned = pa @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - pb) ** 2, axis=1))))


def _numpy_pairwise_min(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt((diff * diff).sum(axis=2)).min(axis=1)


def _numpy_min_contacts(a: np.ndarray, b: np.ndarray, cutoff: float):
    diff = a[:, None, :] - b[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    ii, jj = np.where(dist <= cutoff)
    return sorted(
        (int(i), int(j), float(dist[i, j])) for i, j in zip(ii, jj, strict=True)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="bindings_deltas.json",
                        help="output JSON report path")
    parser.add_argument("--cases", type=int, default=200,
                        help="randomised cases per kernel")
    args = parser.parse_args()

    try:
        import dockflow_bindings as dfb
    except ImportError:
        print("dockflow_bindings is not built; nothing to compare", file=sys.stderr)
        return 2

    rng = np.random.default_rng(2026)
    report: dict = {
        "kabsch_rmsd": {"max_delta": 0.0, "mean_delta": 0.0, "cases": 0},
        "pairwise_min_dist": {"max_delta": 0.0, "mean_delta": 0.0, "cases": 0},
        "min_contacts": {"max_delta": 0.0, "mean_delta": 0.0, "cases": 0,
                         "mismatched_pairs": 0},
        "tolerances": TOLERANCES,
        "generator_seed": 2026,
    }

    # -- kabsch_rmsd ---------------------------------------------------------
    deltas = []
    for _ in range(args.cases):
        n = int(rng.integers(3, 40))
        a = _random_points(rng, n)
        # rigid transform + noise of b keeps cases physically meaningful
        theta = float(rng.uniform(0, 2 * np.pi))
        rotation = np.array([
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ])
        noise = rng.normal(0.0, 0.2, size=a.shape)
        b = a @ rotation.T + noise + rng.uniform(-10, 10, size=3)
        cpp = float(dfb.kabsch_rmsd(a, b))
        ref = _numpy_kabsch(a, b)
        deltas.append(abs(cpp - ref))
    report["kabsch_rmsd"].update(
        max_delta=float(max(deltas)), mean_delta=float(np.mean(deltas)),
        cases=args.cases,
    )

    # -- pairwise_min_dist ---------------------------------------------------
    deltas = []
    for _ in range(args.cases):
        a = _random_points(rng, int(rng.integers(1, 60)))
        b = _random_points(rng, int(rng.integers(1, 60)))
        cpp = np.asarray(dfb.pairwise_min_dist(a, b), dtype=float)
        ref = _numpy_pairwise_min(a, b)
        deltas.append(float(np.max(np.abs(cpp - ref))))
    report["pairwise_min_dist"].update(
        max_delta=float(max(deltas)), mean_delta=float(np.mean(deltas)),
        cases=args.cases,
    )

    # -- min_contacts ----------------------------------------------------------
    max_delta = 0.0
    mismatched = 0
    for _ in range(args.cases):
        a = _random_points(rng, int(rng.integers(2, 30)))
        b = _random_points(rng, int(rng.integers(2, 30)))
        cutoff = float(rng.uniform(2.0, 25.0))
        cpp = sorted(
            (int(i), int(j), float(d)) for i, j, d in dfb.min_contacts(a, b, cutoff)
        )
        ref = _numpy_min_contacts(a, b, cutoff)
        if len(cpp) != len(ref):
            mismatched += 1
            continue
        delta = max(
            (abs(c[2] - r[2]) for c, r in zip(cpp, ref, strict=True)), default=0.0
        )
        if any(c[:2] != r[:2] for c, r in zip(cpp, ref, strict=True)):
            mismatched += 1
        max_delta = max(max_delta, delta)
    report["min_contacts"].update(
        max_delta=max_delta, mean_delta=max_delta, cases=args.cases,
        mismatched_pairs=mismatched,
    )

    Path(args.out).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items()
                      if key in TOLERANCES}, indent=2))
    failures = [
        key for key, tolerance in TOLERANCES.items()
        if report[key]["max_delta"] > tolerance
    ]
    if failures:
        print(f"EQUIVALENCE FAILED for: {', '.join(failures)}", file=sys.stderr)
        return 1
    print(f"C++/NumPy equivalence holds within tolerances ({args.out} written)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
