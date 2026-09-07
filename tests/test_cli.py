"""CLI tests: parser wiring and offline-safe subcommands."""

from __future__ import annotations

from pathlib import Path

import pytest

from dockflow_core.cli import build_parser, main


def test_parser_help(capsys):
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--help"])
    out = capsys.readouterr().out
    for command in ("download", "prep", "gridbox", "dock", "analyze",
                    "visualize", "run", "info", "gui"):
        assert command in out


def test_parser_version(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    # version follows the installed distribution metadata (pyproject.toml);
    # assert the prefix and that it matches dockflow_core.__version__
    import re

    assert re.search(r"dockflow \d+\.\d+\.\d+", out)
    from dockflow_core import __version__

    assert __version__ in out


def test_main_requires_command():
    with pytest.raises(SystemExit):
        main([])


def test_cli_prep_receptor(receptor_pdb_path: Path, tmp_path: Path, capsys):
    exit_code = main([
        "prep", "receptor",
        "--in", str(receptor_pdb_path),
        "--out-dir", str(tmp_path / "prepared"),
        "--engine", "none",
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "receptor PDBQT" in out
    assert (tmp_path / "prepared" / "receptor.pdbqt").is_file()


def test_cli_prep_receptor_missing_input(tmp_path: Path, capsys):
    exit_code = main([
        "prep", "receptor", "--in", str(tmp_path / "nope.pdb"),
        "--out-dir", str(tmp_path), "--engine", "none",
    ])
    assert exit_code == 2
    assert "not found" in capsys.readouterr().err


def test_cli_gridbox_from_structure(receptor_pdb_path: Path, tmp_path: Path,
                                    capsys):
    out = tmp_path / "box.txt"
    exit_code = main([
        "gridbox", "--structure", str(receptor_pdb_path),
        "--resname", "BEN", "--padding", "3", "--out", str(out),
    ])
    assert exit_code == 0
    text = capsys.readouterr().out
    assert "GridBox(" in text
    content = out.read_text(encoding="utf-8")
    assert "center_x" in content and "size_x" in content


def test_cli_gridbox_explicit(tmp_path: Path):
    out = tmp_path / "box.txt"
    exit_code = main([
        "gridbox", "--center", "1,2,3", "--size", "20,22,24", "--out", str(out),
    ])
    assert exit_code == 0
    assert "center_x = 1.000" in out.read_text(encoding="utf-8")


def test_cli_gridbox_requires_input(tmp_path: Path, capsys):
    with pytest.raises(SystemExit):
        main(["gridbox", "--out", str(tmp_path / "box.txt")])


def test_cli_info(capsys):
    exit_code = main(["info"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "environment report" in out
    assert "dockflow-automator" in out
    assert "python" in out


def test_cli_download_requires_args():
    with pytest.raises(SystemExit):
        main(["download"])


def test_cli_download_bad_pdb_id(tmp_path: Path, capsys):
    # invalid id is rejected before any network access; main() maps the
    # DockFlowError to exit code 2.
    exit_code = main(["download", "pdb", "--id", "XXXX", "--out", str(tmp_path)])
    assert exit_code == 2
    assert "invalid PDB id" in capsys.readouterr().err


def test_cli_analyze(receptor_pdb_path: Path, docked_pdbqt_path: Path,
                     tmp_path: Path, capsys):
    # prepare receptor and analyze a fixture docking output
    main([
        "prep", "receptor", "--in", str(receptor_pdb_path),
        "--out-dir", str(tmp_path / "prep"), "--engine", "none",
    ])
    exit_code = main([
        "analyze", "--docking", str(docked_pdbqt_path),
        "--receptor", str(tmp_path / "prep" / "receptor.pdbqt"),
        "--out-dir", str(tmp_path / "analysis"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "interactions.json" in out
    assert (tmp_path / "analysis" / "interactions.json").is_file()


def test_cli_analyze_missing_dir(tmp_path: Path):
    with pytest.raises(SystemExit):
        main(["analyze", "--docking", str(tmp_path / "empty"),
              "--receptor", "x.pdbqt"])


def test_cli_dock_missing_ligand(tmp_path: Path):
    with pytest.raises(SystemExit):
        main([
            "dock", "--receptor", "r.pdbqt",
            "--ligands", str(tmp_path / "no-such-file.pdbqt"),
            "--center", "0,0,0", "--size", "10,10,10",
        ])


def test_cli_dock_no_search_space(tmp_path: Path, receptor_pdb_path: Path,
                                  ligand_pdbqt_path: Path):
    with pytest.raises(SystemExit):
        main([
            "dock", "--receptor", str(receptor_pdb_path),
            "--ligands", str(ligand_pdbqt_path),
        ])


def test_cli_run_config_missing(tmp_path: Path, capsys):
    exit_code = main(["run", "--config", str(tmp_path / "missing.yaml")])
    assert exit_code != 0


def test_download_command_shape():
    parser = build_parser()
    args = parser.parse_args([
        "download", "ligand", "--pubchem", "aspirin", "--out", "raw",
    ])
    assert args.command == "download" and args.kind == "ligand"
    assert args.pubchem == "aspirin"


def test_prep_ligand_requires_source():
    parser = build_parser()
    args = parser.parse_args(["prep", "ligand", "--out-dir", "x"])
    assert args.smiles is None and args.input is None


# ---------------------------------------------------------------------------
# New CLI surface (audit items 23, 27, 31)
# ---------------------------------------------------------------------------
def test_run_dry_run_prints_resolved_config(tmp_path: Path, capsys):
    config = tmp_path / "run.yaml"
    config.write_text(
        "run_id: dry\nworkdir: dryruns\ntarget:\n  pdb_id: 1HVR\n"
        "ligands:\n  - id: l1\n    smiles: CCO\n"
        "gridbox:\n  source: ligand\n  reference_ligand_resname: XK2\n"
        "docking:\n  backend: auto\n  scoring: vina\n  exhaustiveness: 16\n",
        encoding="utf-8",
    )
    code = main(["run", "--config", str(config), "--dry-run"])
    out = capsys.readouterr().out
    assert code == 0
    assert "dry run" in out
    assert "resolved_engine_chain" in out
    assert "dryruns" in out
    # nothing was executed: no run directory beyond config parsing
    assert not (tmp_path / "dryruns" / "dry" / "manifest.json").exists()


def test_run_parser_accepts_stage_and_force():
    parser = build_parser()
    args = parser.parse_args(
        ["run", "--config", "c.yaml", "--stage", "dock", "--stage", "analyze",
         "--force", "--dry-run"]
    )
    assert args.stage == ["dock", "analyze"]
    assert args.force and args.dry_run


def test_run_rejects_unknown_stage():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["run", "--config", "c.yaml", "--stage", "explode"])


def test_enrich_cli(tmp_path: Path, capsys):
    summary = tmp_path / "summary.csv"
    summary.write_text(
        "ligand,pose,affinity_kcal_mol,rmsd_lb,rmsd_ub,crystal_rmsd,"
        "docking_score_efficiency,num_heavy_atoms,runtime_s,backend,error\n"
        "active1,1,-9.0,0,0,, ,30,1,python,\n"
        "active2,1,-8.5,0,0,, ,30,1,python,\n"
        "decoy1,1,-5.0,0,0,, ,30,1,python,\n"
        "decoy2,1,-4.5,0,0,, ,30,1,python,\n",
        encoding="utf-8",
    )
    actives = tmp_path / "actives.txt"
    actives.write_text("active1\nactive2\n", encoding="utf-8")
    code = main(["enrich", "--summary", str(summary), "--actives", str(actives)])
    out = capsys.readouterr().out
    assert code == 0
    assert "roc_auc" in out and "bedroc" in out
    assert "2 actives / 2 decoys" in out
    # comma-separated actives also work
    code = main(["enrich", "--summary", str(summary), "--actives", "active1,active2",
                 "--json"])
    out = capsys.readouterr().out
    assert code == 0 and '"roc_auc": 1.0' in out


def test_enrich_cli_missing_actives(tmp_path: Path):
    summary = tmp_path / "summary.csv"
    summary.write_text(
        "ligand,pose,affinity_kcal_mol\nx,1,-9.0\n", encoding="utf-8"
    )
    code = main(["enrich", "--summary", str(summary), "--actives", "ghost"])
    assert code == 2


def test_log_format_json_flag():
    parser = build_parser()
    args = parser.parse_args(["--log-format", "json", "info"])
    assert args.log_format == "json"


def test_json_log_formatter_emits_events():
    import json
    import logging

    from dockflow_core.utils import JsonLogFormatter

    record = logging.LogRecord(
        name="dockflow.pipeline", level=logging.WARNING, pathname=__file__,
        lineno=1, msg="engine %s failed", args=("openbabel",), exc_info=None,
    )
    record.run_id = "run_1"
    event = json.loads(JsonLogFormatter().format(record))
    assert event["level"] == "WARNING"
    assert event["message"] == "engine openbabel failed"
    assert event["logger"] == "dockflow.pipeline"
    assert event["run_id"] == "run_1"
    assert "timestamp" in event
