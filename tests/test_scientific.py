"""Scientific validation tests (audit item 6, peer item 16).

Unlike the offline unit suite these tests run **real docking** (AutoDock
Vina python bindings) against **real structures** downloaded from RCSB,
and assert scientific outcomes, not just software behaviour:

* 1HVR redocking pose recovery: the co-crystallized ligand XK2 must be
  recovered within 2.0 A of its crystal pose (conventional redocking
  success criterion), using a fully pinned configuration (rdkit receptor
  engine, single-threaded docking, seed 2026, exhaustiveness 8) so the
  assertion is reproducible on any machine;
* seed determinism: the same seed must reproduce the same scores.

They are excluded from the default run (``-m 'not network and not
scientific'``) and executed in CI by the ``reproducibility`` workflow.
Runtime is minutes, not seconds.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytestmark = [pytest.mark.scientific, pytest.mark.network]

for _module, _dist in (("vina", "vina"), ("meeko", "meeko"), ("rdkit", "rdkit")):
    try:
        __import__(_module)
    except ImportError:  # pragma: no cover - scientific stack missing
        pytest.skip(f"{_dist} not installed: scientific tests require the "
                    "full stack (pip install 'dockflow-automator[prep,engine]')",
                    allow_module_level=True)


def _config(workdir: Path, exhaustiveness: int):
    """Pinned, machine-independent validation configuration.

    Every knob that changes the *science* is pinned here so the assertion
    means the same thing on a laptop, a 4-vCPU GitHub runner and a 64-core
    workstation:

    * ``receptor.engine="rdkit"``: the one hydration engine available on
      every platform via the ``prep`` extra (openbabel-wheel is Linux-only,
      so ``auto`` would silently resolve differently per machine - the
      CI failure that motivated this pin).
    * ``docking.cpu=1``: Vina results depend on the number of threads
      (same seed + different cpu counts -> different search trajectories).
      Single-threaded docking is deterministic on every machine.
    """
    from dockflow_core.pipeline import PipelineConfig

    return PipelineConfig(
        workdir=workdir,
        run_id="scientific_1hvr",
        target={"pdb_id": "1HVR"},
        ligands=[{"id": "xk2", "pdb_ligand": "XK2"}],
        receptor={"engine": "rdkit", "charge_model": "gasteiger"},
        gridbox={"source": "ligand", "reference_ligand_resname": "XK2",
                 "padding": 4.0},
        docking={"backend": "auto", "scoring": "vina",
                 "exhaustiveness": exhaustiveness, "num_modes": 9,
                 "seed": 2026, "cpu": 1, "timeout": 1200},
        analysis={"top_poses": 3},
        visualization={"enabled": False},
    )


def test_1hvr_redocking_pose_recovery(tmp_path: Path):
    """The 1HVR/XK2 example is a validation, not just a smoke test.

    Asserts the best docked pose lies within 2.0 A (symmetry-tolerant
    heavy-atom Kabsch RMSD) of the co-crystallized XK2 pose.
    """
    from dockflow_core.pipeline import DockingPipeline

    report = DockingPipeline(_config(tmp_path, exhaustiveness=8)).run()
    assert report.ok, report.error
    results = report.docking["results"]
    assert results and results[0]["poses"]
    crystal = [pose["crystal_rmsd"] for pose in results[0]["poses"]
               if pose.get("crystal_rmsd") is not None]
    assert crystal, "crystal RMSD must be computed for the reference ligand"
    best_rmsd = min(crystal)
    assert best_rmsd <= 2.0, (
        f"redocking failed to recover the crystal pose "
        f"(best RMSD {best_rmsd:.2f} A > 2.0 A)"
    )
    # the score must be in a physically sensible range for XK2/HIV-1 PR
    best_affinity = results[0]["best_affinity"]
    assert -15.0 <= best_affinity <= -6.0, best_affinity


def test_vina_seed_determinism(tmp_path: Path):
    """Same seed + same box + same ligand -> identical scores (item 6).

    Vina is reproducible when the seed AND the cpu/thread count are held
    fixed (the test config pins cpu=1, which additionally makes the result
    independent of the machine's core count).  This test pins that property
    for the versions in the environment fingerprint, and is the local
    counterpart of the CI reproducibility workflow (byte-compare of two
    full runs).
    """
    import csv
    import io

    from dockflow_core.pipeline import DockingPipeline

    def _scientific_rows(text: str) -> list[dict]:
        """summary.csv rows without the wall-clock runtime column."""
        return [
            {key: value for key, value in row.items() if key != "runtime_s"}
            for row in csv.DictReader(io.StringIO(text))
        ]

    runs = []
    for attempt in range(2):
        workdir = tmp_path / f"run{attempt}"
        config = _config(workdir, exhaustiveness=1)
        config.run_id = "scientific_1hvr"
        report = DockingPipeline(config).run()
        assert report.ok, report.error
        summary = (report.run_dir / "docking" / "summary.csv").read_text(
            encoding="utf-8"
        )
        rows = _scientific_rows(summary)
        assert rows and rows[0]["affinity_kcal_mol"]  # rows actually parsed
        runs.append(rows)
    assert runs[0] == runs[1], (
        "same seed produced different scientific columns in summary.csv"
    )
