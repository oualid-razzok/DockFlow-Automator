"""End-to-end pipeline test with mocked downloads and docking engine."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dockflow_core import pipeline as pipeline_mod
from dockflow_core.models import DockingResult, LigandRecord, PoseRecord
from dockflow_core.pipeline import (
    DockingPipeline,
    PipelineConfig,
    PipelineEvents,
)


def _write_pdbqt_fixture(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def config(tmp_path: Path, receptor_pdb_path: Path, ligand_pdbqt_path: Path) -> PipelineConfig:
    return PipelineConfig(
        workdir=tmp_path / "runs",
        run_id="unit_run",
        target={"file": str(receptor_pdb_path)},
        ligands=[{"id": "lig", "pdbqt": str(ligand_pdbqt_path)}],
        receptor={"engine": "none", "charge_model": "zero"},
        gridbox={"source": "residues", "chain": "A", "residues": [1, 2, 3],
                 "padding": 4.0},
        docking={"backend": "auto", "exhaustiveness": 4, "seed": 1,
                 "timeout": 60},
        analysis={"top_poses": 2},
        visualization={"enabled": False},
    )


def _mock_docking(monkeypatch, tmp_path: Path, docked_pdbqt_text: str):
    """Replace VinaEngine with a fake engine writing real output files."""

    class FakeVinaEngine:
        backend = "mocked"

        def __init__(self, config, backend="auto", workdir=None,
                     vina_exec=None, smina_exec=None, gnina_exec=None):
            pass

        def dock_batch(self, receptor, ligand_pdbqts, out_dir=None,
                       ligand_records=None, parallel=1, progress=None,
                       stop_event=None, on_result=None):
            out_dir = Path(out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            results = []
            for index, ligand_path in enumerate(ligand_pdbqts):
                out_path = out_dir / f"{Path(ligand_path).stem}_out.pdbqt"
                out_path.write_text(docked_pdbqt_text, encoding="utf-8")
                record = (ligand_records[index] if ligand_records else None) or \
                    LigandRecord(identifier=Path(ligand_path).stem)
                result = DockingResult(
                    ligand=record,
                    ligand_name=record.identifier,
                    poses=[
                        PoseRecord(model=1, affinity=-9.423),
                        PoseRecord(model=2, affinity=-8.711),
                        PoseRecord(model=3, affinity=-7.905),
                    ],
                    out_path=out_path,
                    log_text="fake",
                    runtime=0.01,
                    backend="mocked",
                )
                results.append(result)
                if on_result is not None:
                    on_result(result)
                if progress:
                    progress((index + 1) / len(ligand_pdbqts), "mocked")
            return results

    monkeypatch.setattr(pipeline_mod, "VinaEngine", FakeVinaEngine)


def test_pipeline_end_to_end(monkeypatch, config, docked_pdbqt_text,
                             ligand_pdbqt_path, tmp_path):
    _mock_docking(monkeypatch, tmp_path, docked_pdbqt_text)
    events = PipelineEvents()
    steps: list[tuple] = []
    events.on_step = lambda step, status, detail: steps.append((step, status))
    logs: list[str] = []
    events.on_log = logs.append
    pipeline = DockingPipeline(config, events=events)
    report = pipeline.run()

    assert report.ok, report.error
    assert report.run_id == "unit_run"
    run_dir = report.run_dir
    assert run_dir is not None and run_dir.is_dir()
    # directory structure
    for sub in ("raw", "prepared", "docking", "analysis", "logs"):
        assert (run_dir / sub).is_dir()
    # artefacts
    assert (run_dir / "manifest.json").is_file()
    assert (run_dir / "report.md").is_file()
    assert (run_dir / "gridbox.txt").is_file()
    assert (run_dir / "prepared" / "receptor.pdbqt").is_file()
    assert (run_dir / "docking" / "summary.csv").is_file()
    assert (run_dir / "analysis" / "interactions.json").is_file()
    # manifest content
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["ok"] is True
    assert manifest["docking"]["backend"] == "mocked"
    assert manifest["docking"]["num_ok"] == 1
    assert manifest["gridbox"]["center"]
    # report content
    report_md = (run_dir / "report.md").read_text(encoding="utf-8")
    assert "# DockFlow-Automator run report" in report_md
    assert "lig" in report_md and "-9.42" in report_md
    # events were emitted
    assert ("download", "done") in steps
    assert ("docking", "done") in steps
    assert ("report", "done") in steps
    assert any("grid box" in entry for entry in logs)


def test_pipeline_config_validation():
    with pytest.raises(Exception, match="target"):
        PipelineConfig.from_dict({"ligands": [{"smiles": "CCO"}]})
    with pytest.raises(Exception, match="ligand"):
        PipelineConfig.from_dict({"target": {"pdb_id": "1HVR"}})
    with pytest.raises(Exception, match="unrecognised ligand entry|recognised"):
        PipelineConfig.from_dict({
            "target": {"pdb_id": "1HVR"},
            "ligands": [{"bogus": 1}],
        })
    with pytest.raises(Exception, match="scoring"):
        PipelineConfig.from_dict({
            "target": {"pdb_id": "1HVR"},
            "ligands": [{"smiles": "CCO"}],
            "docking": {"scoring": "not-a-forcefield"},
        })


def test_pipeline_config_yaml_roundtrip(tmp_path: Path, config):
    import yaml

    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config.to_dict()), encoding="utf-8")
    again = PipelineConfig.from_yaml(path)
    assert again.target == config.target
    assert again.ligands == config.ligands
    assert again.gridbox == config.gridbox


def test_pipeline_cancellation(config, docked_pdbqt_text, ligand_pdbqt_path,
                               tmp_path):
    pipeline = DockingPipeline(config)
    # cancel before the run starts
    pipeline.cancel()
    report = pipeline.run()
    assert report.cancelled or not report.ok


def test_pipeline_reports_failures(config, tmp_path, monkeypatch):
    # a ligand path that does not exist -> pipeline fails cleanly
    from dockflow_core.downloader import DownloadError

    bad = PipelineConfig(
        workdir=tmp_path / "runs",
        run_id="bad_run",
        target={"pdb_id": "1HVR"},
        ligands=[{"smiles": "CCO"}],
        visualization={"enabled": False},
    )

    def fail_download(self, spec, out):
        raise DownloadError("network down (mock)")

    monkeypatch.setattr(pipeline_mod.TargetResolver, "resolve", fail_download)
    report = DockingPipeline(bad).run()
    assert not report.ok
    assert report.error is not None
    assert "network down" in report.error
    assert report.run_dir is not None
    assert (report.run_dir / "manifest.json").is_file()


# ---------------------------------------------------------------------------
# Stage selection, checkpoint/resume and provenance (items 4, 5, 27, 28)
# ---------------------------------------------------------------------------
def test_pipeline_records_stages_and_duration(config, monkeypatch,
                                              docked_pdbqt_text):
    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    pipeline = DockingPipeline(config)
    report = pipeline.run()
    assert report.ok, report.error
    manifest = json.loads((report.run_dir / "manifest.json").read_text(encoding="utf-8"))
    # audit item 4: duration + per-stage audit trail
    assert manifest["duration_s"] is not None and manifest["duration_s"] >= 0.0
    names = [stage["name"] for stage in manifest["stages"]]
    for expected in ("download", "prepare_receptor", "prepare_ligands", "gridbox",
                     "docking", "analysis", "report"):
        assert expected in names
    statuses = {stage["name"]: stage["status"] for stage in manifest["stages"]}
    assert all(status == "done" for name, status in statuses.items()
               if name != "visualization")
    assert statuses.get("visualization") == "skipped"
    for stage in manifest["stages"]:
        assert stage["started_at"] and stage["ended_at"]
        assert stage["duration_s"] is not None
    # audit item 5: environment fingerprint + sidecar file
    assert manifest["environment"]["dockflow"]
    sidecar = report.run_dir / "environment.json"
    assert sidecar.is_file()
    environment = json.loads(sidecar.read_text(encoding="utf-8"))
    assert environment["dockflow"] == manifest["environment"]["dockflow"]


def test_pipeline_manifest_records_receptor_decisions(config, monkeypatch,
                                                       docked_pdbqt_text):
    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    report = DockingPipeline(config).run()
    assert report.ok
    manifest = json.loads((report.run_dir / "manifest.json").read_text(encoding="utf-8"))
    receptor = manifest["receptor"]
    # audit item 10: every silent preparation decision, explicit
    decisions = receptor["decisions"]
    for key in ("rigid_receptor", "protonation", "tautomer", "hydrogen_addition",
                "charge_model", "waters", "metals_cofactors",
                "alternate_conformations", "chains"):
        assert key in decisions, f"missing decision: {key}"
    # audit item 11: removed species recorded (ZN2 + BEN + waters)
    assert receptor["removed"].get("HOH") == 2
    assert receptor["removed"].get("ZN2") == 1


def test_pipeline_single_stage_run_skips_others(config, monkeypatch,
                                                docked_pdbqt_text):
    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    pipeline = DockingPipeline(config)
    report = pipeline.run()
    assert report.ok
    # now re-run only the analysis stage from the same run directory
    calls = {"docking": 0}

    class CountingEngine:
        backend = "mocked"

        def __init__(self, *args, **kwargs):
            pass

        def dock_batch(self, *args, **kwargs):
            calls["docking"] += 1
            raise AssertionError("docking must not re-run for --stage analyze")

    monkeypatch.setattr(pipeline_mod, "VinaEngine", CountingEngine)
    pipeline2 = DockingPipeline(config)
    report2 = pipeline2.run(stages=["analyze"])
    assert report2.ok, report2.error
    assert calls["docking"] == 0
    manifest = json.loads(
        (report2.run_dir / "manifest.json").read_text(encoding="utf-8")
    )
    stages = {stage["name"]: stage["status"] for stage in manifest["stages"]}
    assert stages["docking"] == "skipped"
    assert stages["analysis"] == "done"
    assert (report2.run_dir / "report.md").is_file()  # report always written


def test_pipeline_checkpoint_resume_skips_completed(config, monkeypatch,
                                                     docked_pdbqt_text,
                                                     ligand_pdbqt_path):
    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    report = DockingPipeline(config).run()
    assert report.ok
    progress = json.loads((report.run_dir / "progress.json").read_text(encoding="utf-8"))
    assert "ligand" in progress["completed"]  # keyed by pdbqt stem
    assert progress["completed"]["ligand"]["status"] == "ok"
    assert (report.run_dir / "docking" / "ligand_out.pdbqt").is_file()

    # second invocation: completed ligands must be skipped (no re-docking)
    dock_calls = {"n": 0}

    class NoDockEngine:
        backend = "mocked"

        def __init__(self, *args, **kwargs):
            pass

        def dock_batch(self, receptor, ligand_pdbqts, **kwargs):
            dock_calls["n"] += len(ligand_pdbqts)
            return []

    monkeypatch.setattr(pipeline_mod, "VinaEngine", NoDockEngine)
    config2 = PipelineConfig(**{**config.to_dict()})
    report2 = DockingPipeline(config2).run()
    assert report2.ok, report2.error
    assert dock_calls["n"] == 0  # lig already checkpointed
    manifest = json.loads(
        (report2.run_dir / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["docking"]["resumed"] == 1

    # force re-docks everything (audit item 28: --force)
    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    report3 = DockingPipeline(config2).run(force=True)
    assert report3.ok
    assert json.loads((report3.run_dir / "progress.json").read_text(
        encoding="utf-8"))["completed"]["ligand"]["status"] == "ok"


def test_pipeline_rejects_unknown_stage(config):
    with pytest.raises(Exception, match="unknown stage"):
        DockingPipeline(config).run(stages=["explode"])


def test_pipeline_gridbox_warnings_in_manifest(config, monkeypatch,
                                               docked_pdbqt_text):
    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    # a huge explicit box must be refused (audit item 12)
    huge = PipelineConfig(**{**config.to_dict()})
    huge.gridbox = {"source": "explicit", "center": [0, 0, 0],
                    "size": [40, 40, 40]}
    pipeline = DockingPipeline(huge)
    report = pipeline.run()
    assert not report.ok
    assert "exceeds AutoDock Vina's hard cap" in report.error


# ---------------------------------------------------------------------------
# Integrity checksums (work-order item 29)
# ---------------------------------------------------------------------------
def test_pipeline_manifest_has_input_output_checksums(config, monkeypatch,
                                                       docked_pdbqt_text):
    """SHA-256 of key files, run-dir-relative, verifiable afterwards."""
    import hashlib as _hashlib

    _mock_docking(monkeypatch, config.workdir, docked_pdbqt_text)
    report = DockingPipeline(config).run()
    assert report.ok
    manifest = json.loads((report.run_dir / "manifest.json").read_text(encoding="utf-8"))
    checksums: dict = manifest["checksums"]
    assert checksums, "no checksums recorded"
    # run-dir-relative keys only (portable manifest)
    assert all(not key.startswith("/") for key in checksums)
    # the docking protocol + results must be covered
    assert "gridbox.txt" in checksums
    assert any(key.startswith("prepared/") for key in checksums)
    # (raw/ inputs appear only for downloaded targets; this fixture is local)
    # spot-verify one entry against an independent re-hash
    sample = "gridbox.txt"
    digest = _hashlib.sha256(
        (report.run_dir / sample).read_bytes()).hexdigest()
    assert checksums[sample] == digest
