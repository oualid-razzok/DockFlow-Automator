"""Docking engine tests: config, CLI args, log parsing, ranking, backends."""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from dockflow_core.docker_engine import (
    VinaConfig,
    detect_backends,
    parse_vina_log,
    rank_results,
    write_summary_csv,
)
from dockflow_core.gridbox import GridBox
from dockflow_core.models import DockingResult, LigandRecord, PoseRecord


def test_vina_config_cli_args(tmp_path: Path):
    config = VinaConfig(
        center=(1.5, -2.5, 3.25), size=(22, 24, 20), exhaustiveness=16,
        num_modes=5, seed=42, cpu=8, scoring="vinardo",
    )
    args = config.cli_args(
        Path("r.pdbqt"), Path("l.pdbqt"), Path("out.pdbqt"), Path("run.log")
    )
    joined = " ".join(args)
    assert "--receptor r.pdbqt" in joined
    assert "--center_x 1.5000" in joined
    assert "--center_y -2.5000" in joined
    assert "--size_z 20.0000" in joined
    assert "--exhaustiveness 16" in joined
    assert "--num_modes 5" in joined
    assert "--seed 42" in joined
    assert "--cpu 8" in joined
    assert "--scoring vinardo" in joined
    assert "--score_only" not in joined


def test_vina_config_score_only_args(tmp_path: Path):
    config = VinaConfig(center=(0, 0, 0), size=(10, 10, 10))
    args = config.cli_args(Path("r"), Path("l"), Path("o"), Path("log"),
                           mode="score_only")
    assert "--score_only" in args
    assert "--exhaustiveness" not in args


def test_vina_config_from_gridbox():
    box = GridBox(center=(5, 6, 7), size=(20, 30, 40))
    config = VinaConfig.from_gridbox(box, exhaustiveness=32)
    assert config.center == (5.0, 6.0, 7.0)
    assert config.size == (20.0, 30.0, 40.0)
    assert config.exhaustiveness == 32


def test_vina_config_from_config_file(tmp_path: Path):
    config_path = tmp_path / "box.txt"
    config_path.write_text(
        "center_x = 1.0\ncenter_y = 2.0\ncenter_z = 3.0\n"
        "size_x = 20\nsize_y = 22\nsize_z = 24\n"
        "exhaustiveness = 12\n",
        encoding="utf-8",
    )
    config = VinaConfig.from_vina_config_file(config_path)
    assert config.center == (1.0, 2.0, 3.0)
    assert config.exhaustiveness == 12


def test_parse_vina_log(vina_log_text: str):
    rows = parse_vina_log(vina_log_text)
    assert len(rows) == 3
    assert rows[0]["affinity"] == -9.423
    assert rows[1]["rmsd_lb"] == 1.234 and rows[1]["rmsd_ub"] == 2.100
    assert rows[2]["mode"] == 3


def test_parse_vina_log_empty():
    assert parse_vina_log("nothing here\n") == []


def test_rank_results():
    good = DockingResult(ligand_name="good",
                         poses=[PoseRecord(affinity=-9.0)])
    better = DockingResult(ligand_name="better",
                           poses=[PoseRecord(affinity=-11.0)])
    failed = DockingResult(ligand_name="failed", error="boom")
    ranked = rank_results([good, better, failed])
    assert [r.ligand_name for r in ranked] == ["better", "good", "failed"]


def test_write_summary_csv(tmp_path: Path):
    results = [
        DockingResult(
            ligand_name="lig1",
            poses=[PoseRecord(model=1, affinity=-8.5, rmsd_lb=0.0, rmsd_ub=0.0)],
            runtime=12.5,
            backend="cli",
        ),
        DockingResult(ligand_name="lig2", error="failed to dock"),
    ]
    path = write_summary_csv(results, tmp_path / "summary.csv")
    content = path.read_text(encoding="utf-8")
    assert "ligand,pose,affinity_kcal_mol" in content
    assert "lig1,1,-8.500" in content
    assert "failed to dock" in content


def test_docking_result_properties():
    result = DockingResult(
        ligand=LigandRecord(identifier="x"),
        poses=[PoseRecord(affinity=-7.2), PoseRecord(affinity=-6.1)],
    )
    assert result.ligand_name == "x"
    assert result.ok
    assert result.best_affinity == -7.2
    failed = DockingResult(ligand_name="y", error="x")
    assert not failed.ok and failed.best_affinity is None


# ---------------------------------------------------------------------------
# Backend detection with fake executables
# ---------------------------------------------------------------------------
def test_detect_backends_fake_cli(monkeypatch, tmp_path: Path):
    fake_vina = _make_fake_vina(tmp_path)
    _block_python_vina(monkeypatch)
    reports = detect_backends(str(fake_vina), None)
    by_name = {r.name: r for r in reports}
    assert by_name["cli"].available
    assert "exit 0" in by_name["cli"].version or by_name["cli"].version
    assert not by_name["python"].available


def test_engine_no_backend_raises(monkeypatch, tmp_path: Path):
    from dockflow_core.docker_engine import DockingEngineError, VinaEngine

    monkeypatch.setattr("dockflow_core.docker_engine.which", lambda name: None)
    _block_python_vina(monkeypatch)
    with pytest.raises(DockingEngineError, match="no Vina backend"):
        VinaEngine(VinaConfig(), backend="auto", workdir=tmp_path)


def _make_fake_vina(tmp_path: Path) -> Path:

    fake_vina = tmp_path / "vina"
    fake_vina.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    fake_vina.chmod(fake_vina.stat().st_mode | stat.S_IEXEC)
    return fake_vina


def _block_python_vina(monkeypatch):
    """Make ``import vina`` fail so the CLI backend is selected."""
    import builtins

    real_import = builtins.__import__

    def no_vina_import(name, *args, **kwargs):
        if name == "vina":
            raise ImportError("no vina python bindings in this test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_vina_import)


def test_engine_dock_missing_receptor(monkeypatch, tmp_path: Path):
    """A backend is available (fake vina) but input files are missing."""
    from dockflow_core.docker_engine import DockingEngineError, VinaEngine

    fake_vina = _make_fake_vina(tmp_path)
    _block_python_vina(monkeypatch)
    engine = VinaEngine(VinaConfig(), backend="auto", workdir=tmp_path,
                        vina_exec=str(fake_vina))
    assert engine.backend == "cli"
    with pytest.raises(DockingEngineError, match="not found"):
        engine.dock(tmp_path / "missing.pdbqt", tmp_path / "missing2.pdbqt")


def test_engine_dock_with_mocked_cli(monkeypatch, tmp_path: Path,
                                      docked_pdbqt_text, ligand_pdbqt_path,
                                      receptor_pdb_path):
    """Full CLI-backend flow with a fake vina executable that copies a fixture."""
    fake_vina = _make_fake_vina(tmp_path)

    from dockflow_core import docker_engine as engine_mod
    from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions

    prep = ReceptorPreparator(
        ReceptorPrepOptions(engine="none", charge_model="zero")
    ).prepare(receptor_pdb_path, tmp_path)

    def fake_run_command(command, timeout=None, cwd=None, env=None, on_output=None):
        from dockflow_core.utils import CommandResult

        # the CLI writes the poses file and log per the --out/--log args
        out_index = command.index("--out") + 1
        log_index = command.index("--log") + 1
        Path(command[out_index]).write_text(docked_pdbqt_text, encoding="utf-8")
        Path(command[log_index]).write_text("mode | affinity\n", encoding="utf-8")
        return CommandResult(command=command, returncode=0, stdout="done")

    _block_python_vina(monkeypatch)
    monkeypatch.setattr(engine_mod, "run_command", fake_run_command)

    engine = engine_mod.VinaEngine(
        VinaConfig(center=(12, 9, 8), size=(24, 24, 24)),
        backend="auto", workdir=tmp_path / "dock",
        vina_exec=str(fake_vina),
    )
    assert engine.backend == "cli"
    ligand = LigandRecord(identifier="lig", pdbqt_path=ligand_pdbqt_path)
    result = engine.dock(prep.pdbqt_path, ligand_pdbqt_path,
                         ligand_record=ligand)
    assert result.ok
    assert result.backend == "cli"
    assert result.best_affinity == -9.423
    assert len(result.poses) == 3
    assert result.out_path is not None and result.out_path.is_file()
    assert result.log_path is not None and result.log_path.is_file()


def test_engine_batch_progress_and_cancel(monkeypatch, tmp_path: Path,
                                           docked_pdbqt_text, ligand_pdbqt_path,
                                           receptor_pdb_path):
    fake_vina = _make_fake_vina(tmp_path)

    from dockflow_core import docker_engine as engine_mod
    from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions

    prep = ReceptorPreparator(
        ReceptorPrepOptions(engine="none", charge_model="zero")
    ).prepare(receptor_pdb_path, tmp_path)

    def fake_run_command(command, timeout=None, cwd=None, env=None, on_output=None):
        from dockflow_core.utils import CommandResult

        out_index = command.index("--out") + 1
        Path(command[out_index]).write_text(docked_pdbqt_text, encoding="utf-8")
        log_path = Path(command[out_index]).with_suffix(".log")
        log_path.write_text("mode | affinity\n", encoding="utf-8")
        return CommandResult(command=command, returncode=0)

    _block_python_vina(monkeypatch)
    monkeypatch.setattr(engine_mod, "run_command", fake_run_command)

    engine = engine_mod.VinaEngine(
        VinaConfig(center=(12, 9, 8), size=(24, 24, 24)),
        backend="cli", workdir=tmp_path, vina_exec=str(fake_vina),
    )
    progress_calls = []
    ligands = [tmp_path / f"lig{i}.pdbqt" for i in range(3)]
    for path in ligands:
        path.write_text(ligand_pdbqt_path.read_text(encoding="utf-8"), encoding="utf-8")
    results = engine.dock_batch(
        prep.pdbqt_path, ligands, out_dir=tmp_path / "batch",
        progress=lambda f, m: progress_calls.append((f, m)),
    )
    assert len(results) == 3
    assert all(r.ok for r in results)
    assert progress_calls[-1][0] == 1.0

    import threading

    stop = threading.Event()
    stop.set()
    results = engine.dock_batch(prep.pdbqt_path, ligands,
                                 out_dir=tmp_path / "batch2", stop_event=stop)
    assert all(r.error == "cancelled" for r in results)
