"""GNINA backend tests (audit item 32) using a fake executable."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dockflow_core.docker_engine import VinaConfig, VinaEngine, detect_backends
from dockflow_core.pipeline import PipelineConfig


def _make_fake_gnina(tmp_path: Path) -> Path:
    """A cross-platform fake ``gnina`` executable (see test_engine.py)."""
    if sys.platform == "win32":
        fake = tmp_path / "gnina.bat"
        fake.write_text("@echo off\r\necho GNINA 1.0.0\r\n", encoding="utf-8")
        return fake
    fake = tmp_path / "gnina"
    fake.write_text("#!/bin/sh\necho 'GNINA version 1.0.0'\n", encoding="utf-8")
    fake.chmod(0o755)
    return fake


def test_detect_backends_gnina(tmp_path: Path):
    fake = _make_fake_gnina(tmp_path)
    reports = {r.name: r for r in detect_backends(None, None, str(fake))}
    assert reports["gnina"].available
    assert "gnina" in (reports["gnina"].version or "").lower() or reports["gnina"].version


def test_engine_selects_gnina_backend(tmp_path: Path):
    fake = _make_fake_gnina(tmp_path)
    config = VinaConfig(center=(0, 0, 0), size=(20, 20, 20), scoring="vinardo")
    engine = VinaEngine(config, backend="gnina", workdir=tmp_path, gnina_exec=str(fake))
    assert engine.backend == "gnina"


def test_engine_gnina_requested_but_missing(tmp_path: Path):
    config = VinaConfig()
    with pytest.raises(Exception, match="gnina.*not found"):
        VinaEngine(config, backend="gnina", workdir=tmp_path,
                   gnina_exec=str(tmp_path / "missing_gnina"))


def test_gnina_args_shape():
    config = VinaConfig(
        center=(1.0, 2.0, 3.0), size=(20.0, 20.0, 20.0), scoring="vinardo",
        exhaustiveness=16, num_modes=9, seed=7, cpu=2,
        cnn_scoring="rescore", cnn="dense",
    )
    args = config.gnina_args(Path("rec.pdbqt"), Path("lig.pdbqt"), Path("out.pdbqt"))
    joined = " ".join(args)
    for expected in ("--receptor", "--ligand", "--out", "--center_x", "--size_z",
                     "--exhaustiveness 16", "--num_modes 9", "--seed 7", "--cpu 2",
                     "--scoring vinardo", "--cnn_scoring rescore", "--cnn dense"):
        assert expected in joined, joined
    assert "--log" not in joined  # gnina logs to stdout, not --log
    assert "--refine" not in joined


def test_gnina_args_skip_default_scoring():
    config = VinaConfig(scoring="vina")
    args = config.gnina_args(Path("r"), Path("l"), Path("o"))
    assert "--scoring" not in " ".join(args)


def test_pipeline_config_cnn_validation():
    config = PipelineConfig.from_dict({
        "target": {"pdb_id": "1HVR"},
        "ligands": [{"smiles": "CCO"}],
        "docking": {"backend": "gnina", "cnn_scoring": "rescore"},
    })
    assert config.docking["cnn_scoring"] == "rescore"
    # cnn_scoring without the gnina backend must be rejected
    with pytest.raises(Exception, match="requires docking.backend: gnina"):
        PipelineConfig.from_dict({
            "target": {"pdb_id": "1HVR"},
            "ligands": [{"smiles": "CCO"}],
            "docking": {"backend": "auto", "cnn_scoring": "rescore"},
        })
    with pytest.raises(Exception, match="rescore, all or none"):
        PipelineConfig.from_dict({
            "target": {"pdb_id": "1HVR"},
            "ligands": [{"smiles": "CCO"}],
            "docking": {"backend": "gnina", "cnn_scoring": "bogus"},
        })
