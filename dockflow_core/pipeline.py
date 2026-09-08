"""End-to-end automated docking pipeline.

Inspired by the automated workflow logic of the
``omicscodeathon/anticrcwu`` pipeline scripts, this module chains every
DockFlow stage into one runnable object with progress callbacks and
cancellation::

    from dockflow_core.pipeline import DockingPipeline, PipelineConfig

    pipeline = DockingPipeline(PipelineConfig.from_yaml("run.yaml"))
    report = pipeline.run()

Run directory layout (created automatically)::

    <workdir>/<run_id>/
        manifest.json          # machine-readable run summary
        report.md              # human-readable report
        raw/                   # downloaded target + ligands
        prepared/              # receptor.pdbqt, ligand PDBQTs
        gridbox.txt            # vina config of the search space
        docking/               # *_out.pdbqt + logs + summary.csv
        analysis/              # contacts CSV/JSON, residue hotspots
        visualization/         # rendered PNGs + .pse sessions
        logs/                  # pipeline log file
"""

from __future__ import annotations

import json
import os
import threading
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from shutil import which
from typing import Any

import yaml

from . import __version__
from .analyzer import analyze_docking_result
from .config import AppConfig, get_config
from .docker_engine import VinaConfig, VinaEngine, rank_results, write_summary_csv
from .downloader import (
    LigandDownloader,
    PDBDownloader,
    TargetResolver,
    TargetSpec,
)
from .gridbox import (
    GridBox,
    assumption_strength,
    box_from_pocket,
    box_from_residues,
    box_from_structure,
    box_from_vina_config,
    box_to_vina_config,
    validate_box,
)
from .models import LigandRecord
from .preparator import (
    LigandPreparator,
    LigandPrepOptions,
    ReceptorPreparator,
    ReceptorPrepOptions,
    engine_comparability,
)
from .provenance import StageRecorder, scientific_environment, write_environment_sidecar
from .utils import DockFlowError, get_logger, setup_logging, timestamped_run_id

logger = get_logger("pipeline")

__all__ = [
    "PipelineError",
    "PipelineCancelled",
    "PipelineConfig",
    "PipelineEvents",
    "PipelineReport",
    "DockingPipeline",
    "PIPELINE_STEPS",
]

PIPELINE_STEPS = (
    "download",
    "prepare_receptor",
    "prepare_ligands",
    "gridbox",
    "docking",
    "analysis",
    "visualization",
    "report",
)


class PipelineError(DockFlowError):
    """The automated pipeline failed."""


class PipelineCancelled(PipelineError):
    """The user cancelled the pipeline."""


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass
class PipelineConfig:
    """Full configuration of one automated docking run."""

    workdir: str | Path = "dockflow_runs"
    run_id: str | None = None
    target: dict[str, Any] = field(default_factory=lambda: {"pdb_id": None})
    ligands: list[dict[str, Any]] = field(default_factory=list)
    receptor: dict[str, Any] = field(default_factory=lambda: {
        "chains": None, "keep_water": False, "keep_hetero": False,
        "keep_resnames": [], "altloc": "best", "add_hydrogens": True,
        "merge_nonpolar_h": True, "engine": "auto",
    })
    gridbox: dict[str, Any] = field(default_factory=lambda: {
        "source": "auto",  # auto | ligand | residues | explicit
        "padding": 4.0,
        "reference_ligand_resname": None,
        "chain": None,
        "residues": [],
        "center": None,
        "size": None,
    })
    docking: dict[str, Any] = field(default_factory=lambda: {
        "backend": "auto", "exhaustiveness": 8, "num_modes": 9, "refine": 5,
        "seed": None, "cpu": 0, "scoring": "vina", "parallel": 1, "timeout": 3600,
        "export_sdf": True,
        # gnina-only CNN options (audit item 32)
        "cnn_scoring": None, "cnn": None,
    })
    analysis: dict[str, Any] = field(default_factory=lambda: {
        "top_poses": 3, "contacts_cutoff": 5.0,
    })
    visualization: dict[str, Any] = field(default_factory=lambda: {
        "enabled": True, "engine": "auto", "session": True, "top_poses": 5,
    })

    # -- YAML roundtrip -----------------------------------------------------
    @classmethod
    def from_yaml(cls, path: str | Path) -> PipelineConfig:
        text = Path(path).read_text(encoding="utf-8")
        data = yaml.safe_load(text) or {}
        if not isinstance(data, dict):
            raise PipelineError(f"invalid pipeline YAML root in {path}")
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PipelineConfig:
        valid = set(cls.__dataclass_fields__)
        filtered = {k: v for k, v in data.items() if k in valid}
        config = cls(**filtered)
        config._validate()
        return config

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for name in self.__dataclass_fields__:
            value = getattr(self, name)
            if isinstance(value, Path):
                value = str(value)
            result[name] = value
        return result

    def _validate(self) -> None:
        target = self.target or {}
        if not (target.get("pdb_id") or target.get("uniprot") or target.get("file")):
            raise PipelineError(
                "pipeline config needs a target: pdb_id, uniprot or file"
            )
        if not self.ligands:
            raise PipelineError("pipeline config needs at least one ligand entry")
        for entry in self.ligands:
            if not isinstance(entry, dict):
                raise PipelineError(f"ligand entries must be mappings, got {entry!r}")
            if not any(
                entry.get(k) for k in ("smiles", "file", "pubchem", "zinc", "pdb_ligand",
                                       "pdbqt")
            ):
                raise PipelineError(f"ligand entry without a recognised source: {entry!r}")
        scoring = (self.docking or {}).get("scoring", "vina")
        if scoring not in ("vina", "vinardo", "ad4"):
            raise PipelineError("docking.scoring must be vina, vinardo or ad4")
        cnn_scoring = (self.docking or {}).get("cnn_scoring")
        if cnn_scoring is not None:
            if (self.docking or {}).get("backend") not in ("gnina",):
                raise PipelineError(
                    "docking.cnn_scoring requires docking.backend: gnina"
                )
            if cnn_scoring not in ("rescore", "all", "none"):
                raise PipelineError(
                    "docking.cnn_scoring must be rescore, all or none"
                )


# ---------------------------------------------------------------------------
# Events / report
# ---------------------------------------------------------------------------
@dataclass
class PipelineEvents:
    """Callback bundle; every field may be replaced by a callable."""

    on_step: Callable[[str, str, str | None], None] = lambda step, status, detail: None
    on_progress: Callable[[float, str], None] = lambda fraction, message: None
    on_log: Callable[[str], None] = lambda message: None


@dataclass
class PipelineReport:
    """Everything worth knowing after a run."""

    run_id: str = ""
    run_dir: Path | None = None
    ok: bool = False
    error: str | None = None
    cancelled: bool = False
    config: dict[str, Any] = field(default_factory=dict)
    timings: dict[str, float] = field(default_factory=dict)
    # Total wall time of the run in seconds (audit item 4).
    duration_s: float | None = None
    # Per-stage audit trail: wall time, ISO start/stop, exit status,
    # resolved engine/backend (audit item 4).
    stages: list[dict[str, Any]] = field(default_factory=list)
    # Scientific-stack version fingerprint (audit item 5).
    environment: dict[str, Any] = field(default_factory=dict)
    target: dict[str, Any] = field(default_factory=dict)
    receptor: dict[str, Any] = field(default_factory=dict)
    ligands: list[dict[str, Any]] = field(default_factory=list)
    gridbox: dict[str, Any] = field(default_factory=dict)
    docking: dict[str, Any] = field(default_factory=dict)
    analysis: dict[str, Any] = field(default_factory=dict)
    visualization: dict[str, Any] = field(default_factory=dict)
    paths: dict[str, Any] = field(default_factory=dict)
    # SHA-256 of key input/output files, keyed by run-dir-relative path
    # (work-order item 29: tamper detection for published runs).
    checksums: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {}
        for name in self.__dataclass_fields__:
            value = getattr(self, name)
            if isinstance(value, Path):
                value = str(value)
            data[name] = value
        return data

    def save(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return target


# ---------------------------------------------------------------------------
# Degraded-mode banner (final pass, item 8)
# ---------------------------------------------------------------------------
def degraded_mode_banner(decisions: dict, engine: str | None) -> str | None:
    """The unmistakable report banner for scientifically degraded runs.

    Fires when the resolved receptor engine is the dependency-free
    ``none`` fallback (no hydrogens, zero charges, no aromaticity
    perception) or any engine whose results are comparable to nothing.
    Returns the banner line (with the warning sign) or ``None``.
    """
    resolved_engine = (engine or "none").split()[0].strip().lower()
    degraded = (bool(decisions.get("degraded_mode"))
                or resolved_engine == "none"
                or not (decisions.get("comparable_to") or []))
    if not degraded:
        return None
    return (
        "> ⚠ SCIENTIFICALLY DEGRADED MODE: receptor prepared "
        "with the `none` engine (no hydrogens, zero charges, no "
        "aromaticity perception). Results are geometry-only and "
        "NOT recommended for production. Install OpenBabel or "
        "RDKit (`pip install dockflow-automator[prep]`) and re-run."
    )


# ---------------------------------------------------------------------------
# The pipeline
# ---------------------------------------------------------------------------
class DockingPipeline:
    """Orchestrates download -> prepare -> grid box -> dock -> analyze -> render."""

    def __init__(
        self,
        config: PipelineConfig,
        events: PipelineEvents | None = None,
        stop_event: threading.Event | None = None,
        app_config: AppConfig | None = None,
    ) -> None:
        self.config = config
        self.events = events or PipelineEvents()
        self.stop_event = stop_event or threading.Event()
        self.app_config = app_config or get_config()
        self._report = PipelineReport(config=config.to_dict())
        self._run_dir: Path | None = None
        self._stages = StageRecorder()
        # Receptor preparation outcome (kept for the analysis stage: the
        # co-crystallized ligand is extracted from the *filtered* structure
        # for redocking validation).
        self._receptor_result = None

    # -- control ------------------------------------------------------------
    def cancel(self) -> None:
        self.stop_event.set()

    def _check_cancel(self) -> None:
        if self.stop_event.is_set():
            raise PipelineCancelled("cancelled by user")

    def _step(self, name: str, status: str, detail: str | None = None) -> None:
        logger.info("pipeline step %s: %s %s", name, status, detail or "")
        try:
            self.events.on_step(name, status, detail)
        except Exception:  # noqa: BLE001
            logger.debug("on_step callback failed", exc_info=True)

    def _progress(self, fraction: float, message: str) -> None:
        try:
            self.events.on_progress(max(0.0, min(1.0, fraction)), message)
        except Exception:  # noqa: BLE001
            logger.debug("on_progress callback failed", exc_info=True)

    def _log(self, message: str) -> None:
        try:
            self.events.on_log(message)
        except Exception:  # noqa: BLE001
            logger.debug("on_log callback failed", exc_info=True)

    # -- main ---------------------------------------------------------------
    _STAGE_ALIASES = {
        "download": {"download"},
        "prep": {"prepare_receptor", "prepare_ligands"},
        "gridbox": {"gridbox"},
        "dock": {"docking"},
        "analyze": {"analysis"},
        "visualize": {"visualization"},
    }

    def run(
        self,
        stages: Sequence[str] | None = None,
        force: bool = False,
    ) -> PipelineReport:
        """Execute the workflow and return the report.

        Args:
            stages: run only these stages (``download``/``prep``/``gridbox``/
                ``dock``/``analyze``/``visualize``); missing prerequisites are
                loaded from the existing run directory so a single stage can
                be re-run for debugging without redoing the whole pipeline.
            force: ignore the ``progress.json`` checkpoint and re-dock every
                ligand (audit item 28).
        """
        report = self._report
        self._force = force
        requested: set[str] | None = None
        if stages:
            requested = set()
            for stage in stages:
                if stage not in self._STAGE_ALIASES:
                    raise PipelineError(
                        f"unknown stage {stage!r} (valid: "
                        f"{', '.join(sorted(self._STAGE_ALIASES))})"
                    )
                requested |= self._STAGE_ALIASES[stage]
        self._requested_stages = requested
        started = time.perf_counter()
        run_id = self.config.run_id or timestamped_run_id("run")
        report.run_id = run_id
        base = Path(self.config.workdir)
        self._run_dir = base / run_id
        resuming = requested is not None and self._run_dir.is_dir()
        for sub in ("raw", "prepared", "docking", "analysis", "visualization", "logs"):
            (self._run_dir / sub).mkdir(parents=True, exist_ok=True)
        report.run_dir = self._run_dir
        logfile = self._run_dir / "logs" / "pipeline.log"
        setup_logging("INFO", logfile)
        logger.info("DockFlow-Automator %s starting run %s", __version__, run_id)
        self._log(f"run {run_id} -> {self._run_dir}")
        if resuming:
            self._log(f"single-stage run: only {sorted(requested or set())} "
                      "(artifacts loaded from the existing run directory)")
            self._restore_previous_report_state()

        try:
            target_record = self._maybe_run("download", self._download)
            if target_record is None:
                target_record = self._load_previous_target()
            receptor_result = self._maybe_run(
                "prepare_receptor", self._prepare_receptor, target_record,
                on_done=lambda record: setattr(
                    record, "engine",
                    (self._report.receptor or {}).get("engine", "")),
            )
            if receptor_result is None:
                receptor_result = self._load_previous_receptor()
            self._receptor_result = receptor_result
            ligand_records = self._maybe_run(
                "prepare_ligands", self._prepare_ligands,
                on_done=lambda record: setattr(record, "engine", "meeko"),
            )
            if ligand_records is None:
                ligand_records = self._load_previous_ligands()
            box = self._maybe_run("gridbox", self._define_gridbox, target_record,
                                  ligand_records)
            if box is None:
                box = self._load_previous_box()
            results = self._maybe_run(
                "docking", self._dock, receptor_result, ligand_records, box,
                on_done=lambda record: setattr(
                    record, "backend",
                    (self._report.docking or {}).get("backend", "")),
            )
            if results is None:
                results = self._load_previous_results(ligand_records)
            self._maybe_run("analysis", self._analyze, results, receptor_result,
                            target_record, box)
            self._maybe_run("visualization", self._visualize, results,
                            receptor_result, box)
            # The report is the run summary: always written, even for
            # single-stage runs (it is not a --stage-selectable stage).
            with self._stages.stage("report"):
                self._write_report(results, receptor_result, box, target_record)
            report.ok = True
            report.error = None
        except PipelineCancelled as exc:
            report.cancelled = True
            report.error = str(exc)
            report.ok = False
            logger.warning("pipeline cancelled")
        except DockFlowError as exc:
            report.error = str(exc)
            report.ok = False
            logger.error("pipeline failed: %s", exc)
        except Exception as exc:  # noqa: BLE001 - never leak raw tracebacks
            report.error = f"{type(exc).__name__}: {exc}"
            report.ok = False
            logger.exception("unexpected pipeline failure")
        report.timings["total"] = round(time.perf_counter() - started, 2)
        # Provenance (audit items 4 + 5): stage audit trail, total duration,
        # scientific environment fingerprint + sidecar file.
        report.duration_s = report.timings["total"]
        report.stages = self._stages.to_list()
        report.environment = scientific_environment()
        # Input/output integrity (work-order item 29): SHA-256 of the run's
        # key files.  Signing is deliberately NOT implemented - an
        # open-source tool has no key management; hashes give tamper
        # DETECTION (re-hash with `sha256sum -c`-style verification).
        report.checksums = self._checksums()
        manifest = self._run_dir / "manifest.json"
        report.save(manifest)
        sidecar = write_environment_sidecar(self._run_dir, report.environment)
        report.paths["manifest"] = str(manifest)
        report.paths["environment"] = str(sidecar)
        logger.info("run %s finished in %.1fs (ok=%s)", run_id,
                    report.duration_s, report.ok)
        return report

    def _checksums(self) -> dict[str, str]:
        """SHA-256 map of the run's key files (inputs + outputs).

        Keys are run-dir-relative paths so the manifest stays portable;
        missing files are silently skipped (partial runs still publish
        whatever they produced).
        """
        from .utils import sha256_file

        if self._run_dir is None:
            return {}
        run_dir = Path(self._run_dir)
        wanted: list[Path] = []
        # inputs
        wanted += sorted((run_dir / "raw").glob("*.pdb"))
        wanted += sorted((run_dir / "raw").glob("*.sdf"))
        # preparation outputs (the scientific boundary: what Vina consumed)
        wanted += [run_dir / "prepared" / "receptor.pdbqt"]
        wanted += sorted((run_dir / "prepared").glob("*.pdbqt"))
        # protocol + results
        wanted += [run_dir / "gridbox.txt",
                   run_dir / "docking" / "summary.csv"]
        wanted += sorted((run_dir / "docking").glob("*_out.pdbqt"))
        checksums: dict[str, str] = {}
        for path in wanted:
            if path.is_file():
                try:
                    rel = path.relative_to(run_dir).as_posix()
                    checksums[rel] = sha256_file(path)
                except OSError:
                    logger.debug("checksum failed for %s", path, exc_info=True)
        return checksums

    # -- stage selection (audit item 27) ---------------------------------------
    def _maybe_run(self, stage: str, method: Callable[..., Any], *args,
                   on_done: Callable[[Any], None] | None = None, **kwargs):
        """Run ``method`` inside the stage recorder, honouring stage selection.

        Returns ``None`` when the stage was not requested so the caller can
        fall back to loading the previous run's artifacts.  ``on_done``
        receives the :class:`StageRecord` before the stage closes so
        callers can attach the resolved engine/backend (audit item 4).
        """
        requested = getattr(self, "_requested_stages", None)
        if requested is not None and stage not in requested:
            with self._stages.stage(stage) as record:
                record.status = "skipped"
                record.detail = "not requested (--stage)"
            self._step(stage, "skipped", "not requested")
            return None
        with self._stages.stage(stage) as record:
            result = method(*args, **kwargs)
            if on_done is not None:
                try:
                    on_done(record)
                except Exception:  # noqa: BLE001 - provenance must not fail runs
                    logger.debug("on_done stage hook failed", exc_info=True)
            return result

    def _restore_previous_report_state(self) -> None:
        """Restore provenance sections from an existing run's manifest.

        Single-stage re-runs (``--stage analyze``) reload artifacts from
        disk, but the manifest's provenance sections (receptor decisions,
        grid box, ligand prep, docking backend/decisions) are only
        written by the stage that produced them - previously a resumed
        run re-saved a manifest WITHOUT them, silently losing provenance.
        Restoring the sections of stages that are NOT being re-run keeps
        re-analysed runs fully traceable.  Sections of stages that DO
        re-run are overwritten fresh by those stages.
        """
        manifest = self._run_dir / "manifest.json"
        if not manifest.is_file():
            return
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            logger.debug("could not restore previous manifest state",
                         exc_info=True)
            return
        for section in ("target", "receptor", "gridbox", "ligands"):
            if data.get(section) and not getattr(self._report, section, None):
                setattr(self._report, section, data[section])
        docking = data.get("docking") or {}
        for key in ("backend", "num_ligands", "num_ok", "resumed",
                    "warnings", "decisions", "summary_csv"):
            if key in docking and key not in self._report.docking:
                self._report.docking[key] = docking[key]
        for key, value in (data.get("timings") or {}).items():
            self._report.timings.setdefault(key, value)
        if data.get("duration_s") is not None:
            self._report.timings.setdefault("original_run_duration_s",
                                            data["duration_s"])
        logger.info("restored provenance sections from %s", manifest)

    def _load_previous_target(self):
        """Reload the target record from a previous run's artifacts."""
        from .downloader import ProteinRecord

        raw = sorted((self._run_dir / "raw").glob("*.pdb")) if self._run_dir else []
        source_file = self.config.target.get("file")
        if raw:
            record = ProteinRecord(identifier=raw[0].stem, source="pdb (resumed)",
                                   path=raw[0])
        elif source_file and Path(source_file).is_file():
            # file-based target: the input file itself is the artifact
            record = ProteinRecord(identifier=Path(source_file).stem,
                                   source="file (resumed)", path=Path(source_file))
        else:
            raise PipelineError(
                "--stage selection needs an existing run directory with "
                "downloaded artifacts; run the download stage first"
            )
        self._report.target = record.to_dict()
        self._report.paths["target_structure"] = str(record.path)
        return record

    def _load_previous_receptor(self):
        """Reload the receptor preparation result from a previous run."""
        from .preparator import ReceptorPrepResult

        pdbqt = self._run_dir / "prepared" / "receptor.pdbqt"
        if not pdbqt.is_file():
            raise PipelineError(
                "receptor.pdbqt not found in the run directory; run the prep "
                "stage first (dockflow run --stage prep ...)"
            )
        result = ReceptorPrepResult(pdbqt_path=pdbqt, engine="(resumed)")
        if self._report.receptor.get("engine"):
            result.engine = self._report.receptor["engine"]
        self._log(f"receptor loaded from {pdbqt}")
        return result

    def _load_previous_ligands(self) -> list[LigandRecord]:
        """Reload prepared ligands from a previous run."""
        records: list[LigandRecord] = []
        prepared = self._run_dir / "prepared"
        # heavy-atom counts and prep notes are restored from the manifest so
        # efficiency metrics and provenance survive a resumed run
        previous: dict[str, dict] = {}
        manifest = self._run_dir / "manifest.json"
        if manifest.is_file():
            try:
                data = json.loads(manifest.read_text(encoding="utf-8"))
                previous = {
                    entry.get("identifier", ""): entry
                    for entry in data.get("ligands", [])
                    if isinstance(entry, dict)
                }
            except (ValueError, OSError):
                previous = {}
        for index, entry in enumerate(self.config.ligands):
            identifier = entry.get("id") or entry.get("name") or f"ligand_{index + 1}"
            candidates = [prepared / f"{identifier}.pdbqt"]
            if entry.get("pdbqt"):
                candidates.append(Path(entry["pdbqt"]))
            pdbqt = next((path for path in candidates if path.is_file()), None)
            if pdbqt is not None:
                prior = previous.get(identifier, {})
                records.append(LigandRecord(
                    identifier=identifier, source="file",
                    value=str(pdbqt), path=pdbqt, pdbqt_path=pdbqt, status="prepared",
                    num_heavy_atoms=prior.get("num_heavy_atoms"),
                    num_rotatable_bonds=prior.get("num_rotatable_bonds"),
                    prep=prior.get("prep", {}),
                ))
        if not records:
            raise PipelineError(
                "no prepared ligand PDBQTs found; run the prep stage first"
            )
        self._log(f"{len(records)} ligand(s) loaded from {prepared}")
        return records

    def _load_previous_box(self) -> GridBox:
        """Reload the grid box from a previous run's gridbox.txt."""
        config_path = self._run_dir / "gridbox.txt"
        if not config_path.is_file():
            raise PipelineError("gridbox.txt not found; run the gridbox stage first")
        box = box_from_vina_config(config_path)
        # Restore the original derivation (e.g. "pocket:XK2") from the
        # previous manifest so redocking validation keeps working on resume.
        manifest = self._run_dir / "manifest.json"
        if manifest.is_file():
            try:
                data = json.loads(manifest.read_text(encoding="utf-8"))
                source = (data.get("gridbox") or {}).get("source")
                if source:
                    box.source = source
                    box.padding = float((data.get("gridbox") or {}).get("padding", 0.0))
            except (ValueError, OSError):
                logger.debug("could not restore grid box source from manifest")
        base = box.to_dict()
        # Provenance keys that only the original gridbox stage recorded
        # (assumption_strength, warnings, receptor_coverage_pct) survive a
        # resume: they were restored into _report.gridbox by
        # _restore_previous_report_state and are NOT derivable from
        # gridbox.txt.
        extras = {key: value for key, value in (self._report.gridbox or {}).items()
                  if key not in base and key != "vina_config"}
        self._report.gridbox = {**base, **extras,
                                "vina_config": str(config_path)}
        return box

    def _load_previous_results(self, ligand_records):
        """Reload docking results from a previous run's PDBQT outputs."""
        from .docker_engine import detect_backends
        from .models import DockingResult
        from .pdbio import parse_pdbqt_results

        results = []
        for record in ligand_records:
            if record.status != "prepared" or not record.pdbqt_path:
                continue
            out_path = self._run_dir / "docking" / f"{record.pdbqt_path.stem}_out.pdbqt"
            if not out_path.is_file():
                continue
            backend = next(
                (b for b in detect_backends() if b.available), None
            )
            results.append(DockingResult(
                ligand=record, ligand_name=record.identifier,
                poses=parse_pdbqt_results(out_path), out_path=out_path,
                backend=backend.name if backend else "unknown",
            ))
        if not results:
            raise PipelineError(
                "no docking outputs found; run the dock stage first"
            )
        return results

    # ------------------------------------------------------------------ stages
    def _download(self):

        self._check_cancel()
        self._step("download", "running")
        start = time.perf_counter()
        raw_dir = self._run_dir / "raw"
        spec = TargetSpec.parse(self.config.target)
        resolver = TargetResolver(PDBDownloader(cache_dir=self.app_config.cache_dir))
        record = resolver.resolve(spec, raw_dir)
        self._report.target = record.to_dict()
        self._report.paths["target_structure"] = str(record.path)
        self._log(f"target {record.identifier} ({record.source}) -> {record.path}")
        self._step("download", "done", record.identifier)
        self._report.timings["download"] = round(time.perf_counter() - start, 2)
        self._progress(0.08, "target downloaded")
        return record

    def _prepare_receptor(self, target_record):
        self._check_cancel()
        self._step("prepare_receptor", "running")
        start = time.perf_counter()
        options = ReceptorPrepOptions(**{
            key: value
            for key, value in (self.config.receptor or {}).items()
            if key in ReceptorPrepOptions.__dataclass_fields__
        })
        assert target_record.path is not None
        preparator = ReceptorPreparator(options)
        result = preparator.prepare(
            target_record.path, self._run_dir / "prepared", basename="receptor"
        )
        # Flexible side chains (peer item 8): split rigid/flex PDBQTs.
        flex_selection = (self.config.receptor or {}).get("flexible_residues")
        if flex_selection:
            from .flexible import auto_flexible_residues, build_flexible_pdbqts, parse_flexible_residues
            selection = parse_flexible_residues(flex_selection)
            auto_mode = bool(selection) and selection[0].get("auto")
            if auto_mode:
                reference_resname = (self.config.gridbox or {}).get(
                    "reference_ligand_resname") \
                    or (target_record.ligand_codes or [None])[0]
                if reference_resname is None:
                    raise PipelineError(
                        "receptor.flexible_residues: 'auto' needs a reference "
                        "ligand (gridbox.reference_ligand_resname or a "
                        "co-crystallized ligand in the target)")
                selection = auto_flexible_residues(
                    result.pdbqt_path, target_record.path)
                self._log(
                    f"flexible residues (auto): {len(selection)} residues "
                    f"contacting {reference_resname} "
                    "(heuristic: any contact within 4.5 A)")
            if selection:
                try:
                    rigid, flex, report = build_flexible_pdbqts(
                        result.pdbqt_path, selection,
                        self._run_dir / "prepared", basename="receptor")
                    result.rigid_pdbqt_path = rigid
                    result.flex_pdbqt_path = flex
                    result.flexible_residues = report
                    skipped = [r for r in report if r.get("status") != "flexible"]
                    for entry in skipped:
                        result.warnings.append(
                            f"flexible residue {entry['resname']}-"
                            f"{entry['chain']}{entry['resseq']} skipped: "
                            f"{entry['reason']}")
                    result.warnings.append(
                        "flexible-receptor docking is enabled: Vina scores "
                        "are NOT comparable with rigid-receptor scores of "
                        "the same ligand, and the search is slower.")
                except DockFlowError as exc:
                    raise PipelineError(f"flexible residues: {exc}") from exc
        removed_non_water = {
            name: count for name, count in result.removed_resnames.items()
            if name not in {"HOH", "DOD", "WAT", "H2O", "TIP", "TIP3", "SOL"}
        }
        # Comparability of the resolved engine (peer item 4): which other
        # engines would have produced scientifically comparable receptors.
        comparability = engine_comparability(result.engine)
        resolved_base = (result.engine or "none").split()[0].strip().lower()
        degraded = resolved_base == "none" or not comparability["comparable_to"]
        # Every silent preparation decision, made explicit (audit item 10).
        decisions = {
            "rigid_receptor": not bool(result.flex_pdbqt_path),
            "protonation": (
                f"as deposited (PDB); polar hydrogens added by toolkit "
                f"template ({result.engine.split()[0]} engine, pH 7.4 "
                "nominal; NOT a pKa calculation)"
                if options.add_hydrogens else
                "as deposited (no hydrogens added)"
            ),
            "tautomer": "n/a (protein receptor)",
            "hydrogen_addition": (
                f"{result.hydrogens_added} hydrogens added by {result.engine}; "
                + ("non-polar hydrogens merged into heavy atoms"
                   if options.merge_nonpolar_h else "all hydrogens kept explicit")
            ),
            "charge_model": options.charge_model,
            "waters": "kept" if options.keep_water else
                      f"removed ({result.waters_removed} water oxygen atoms)",
            "metals_cofactors": {
                "removed": removed_non_water or "none",
                "kept": options.keep_resnames or
                        ("all hetero" if options.keep_hetero else "none"),
            },
            "alternate_conformations": f"altloc policy: {options.altloc}",
            "chains": options.chains or "all",
            # peer item 4: results must never silently look comparable to
            # runs prepared with a different engine
            "comparable_to": comparability["comparable_to"],
            "incomparable_to": comparability["incomparable_to"],
            # peer items 6, 7, 10: structure quality findings
            "missing_residues": result.missing_residues or "none",
            "disulfides": result.disulfides or "none",
            "reduced_cysteines": result.reduced_cysteines or "none",
            "metal_coordination": result.metal_coordination or "none",
            # peer item 8: flexible side chains
            "flexible_residues": result.flexible_residues or "none",
            # final pass, item 8: machine-readable degraded-mode flag so
            # downstream tools can detect a scientifically degraded run
            # without parsing report prose
            "degraded_mode": degraded,
            "degraded_reason": (
                "receptor prepared with the 'none' engine: no hydrogens "
                "added, zero partial charges, no aromaticity perception - "
                "geometry-only results, not recommended for production"
                if degraded else None
            ),
        }
        self._report.receptor = {
            "pdbqt": str(result.pdbqt_path),
            "pdb": str(result.pdb_path),
            "engine": result.engine,
            "atoms_in": result.atoms_in,
            "atoms_out": result.atoms_out,
            "hydrogens_added": result.hydrogens_added,
            "removed": result.removed_resnames,
            "decisions": decisions,
            "warnings": result.warnings,
        }
        if result.flex_pdbqt_path:
            self._report.receptor.update({
                "rigid_pdbqt": str(result.rigid_pdbqt_path),
                "flex_pdbqt": str(result.flex_pdbqt_path),
            })
        self._report.paths["receptor_pdbqt"] = str(result.pdbqt_path)
        self._log(f"receptor prepared with {result.engine} engine "
                  f"({result.atoms_out} atoms)")
        for warning in result.warnings:
            self._log(f"warning: {warning}")
        self._step("prepare_receptor", "done", f"{result.atoms_out} atoms")
        self._report.timings["prepare_receptor"] = round(time.perf_counter() - start, 2)
        self._progress(0.2, "receptor prepared")
        return result

    def _prepare_ligands(self) -> list[LigandRecord]:
        self._check_cancel()
        self._step("prepare_ligands", "running")
        start = time.perf_counter()
        raw_dir = self._run_dir / "raw"
        prep_dir = self._run_dir / "prepared"
        downloader = LigandDownloader()
        records: list[LigandRecord] = []
        for index, entry in enumerate(self.config.ligands):
            self._check_cancel()
            record = self._resolve_ligand(entry, index, downloader, raw_dir)
            if record.status == "error":
                continue
            if entry.get("pdbqt"):
                # User-supplied pre-prepared ligand.
                record.pdbqt_path = record.path
                record.status = "prepared"
                records.append(record)
                self._log(f"{record.identifier}: using provided PDBQT")
                continue
            try:
                assert record.path is not None
                # Per-ligand preparation options (peer items 11-13):
                # protonation / tautomers / stereochemistry are scientific
                # decisions configured per ligand entry in the YAML.
                options = LigandPrepOptions(**{
                    key: value
                    for key, value in entry.items()
                    if key in LigandPrepOptions.__dataclass_fields__
                })
                preparator = LigandPreparator(options)
                prep = preparator.prepare(record.path, prep_dir, record.identifier)
                record.num_rotatable_bonds = prep.num_rotatable_bonds
                record.num_heavy_atoms = prep.num_heavy_atoms
                # Per-ligand preparation provenance (audit items 10 + 26):
                # hydrogen source, protonation choice, charge model, and
                # the state chemistry of items 11-14.
                prep_dict = {
                    "engine": prep.engine,
                    "charge_model": prep.charge_model,
                    "input_had_explicit_hydrogens": prep.input_had_explicit_hydrogens,
                    "hydrogen_source": (
                        "explicit in input" if prep.input_had_explicit_hydrogens
                        else "added by RDKit/Meeko during preparation"
                    ) if prep.input_had_explicit_hydrogens is not None else "unknown",
                    "protonation": prep.protonation_source,
                    "num_protonation_variants": prep.num_protonation_variants,
                    "tautomer": (prep.tautomers or {}).get(
                        "source", "single input tautomer (no enumeration)"),
                    "warnings": prep.warnings,
                }
                # Full provenance blocks for the manifest (items 11-14):
                if prep.protonation:
                    prep_dict["protonation_detail"] = prep.protonation
                if prep.tautomers:
                    prep_dict["tautomers"] = prep.tautomers
                if prep.stereocenters:
                    prep_dict["stereocenters"] = prep.stereocenters
                if prep.salt_stripping:
                    prep_dict["salt_stripping"] = prep.salt_stripping
                # One LigandRecord per prepared state (peer item 11c: one
                # input ligand -> N PDBQTs -> N summary.csv entries with the
                # state in the ligand id).
                variants = prep.variants or [{
                    "name": record.identifier,
                    "kind": "input", "tag": "", "index": 1,
                    "smiles": prep.smiles, "pdbqt": str(prep.pdbqt_path),
                }]
                for position, variant in enumerate(variants):
                    if position == 0:
                        variant_record = record
                        variant_record.pdbqt_path = Path(variant["pdbqt"])
                    else:
                        variant_record = replace(
                            record, identifier=variant["name"],
                            pdbqt_path=Path(variant["pdbqt"]))
                    variant_record.status = "prepared"
                    variant_record.prep = {
                        **prep_dict, "variant": {
                            key: value for key, value in variant.items()
                            if key != "pdbqt"}}
                    records.append(variant_record)
                self._log(
                    f"{record.identifier}: prepared "
                    f"({len(variants)} state(s), {prep.num_atoms} atoms)")
                for warning in prep.warnings:
                    self._log(f"warning: {record.identifier}: {warning}")
            except DockFlowError as exc:
                record.status = "error"
                record.error = str(exc)
                records.append(record)
                self._log(f"{record.identifier}: preparation failed ({exc})")
            self._progress(
                0.2 + 0.15 * (index + 1) / max(1, len(self.config.ligands)),
                f"ligand {index + 1}/{len(self.config.ligands)}",
            )
        prepared = [r for r in records if r.status == "prepared"]
        if not prepared:
            raise PipelineError("no ligand could be prepared")
        self._report.ligands = [r.to_dict() for r in records]
        self._step("prepare_ligands", "done", f"{len(prepared)} ligand states")
        self._report.timings["prepare_ligands"] = round(time.perf_counter() - start, 2)
        return records

    def _resolve_ligand(
        self, entry: dict[str, Any], index: int, downloader: LigandDownloader, raw_dir: Path
    ) -> LigandRecord:
        identifier = entry.get("id") or entry.get("name") or f"ligand_{index + 1}"
        try:
            if entry.get("smiles"):
                path = raw_dir / f"{identifier}.smi"
                path.write_text(f"{entry['smiles']} {identifier}\n", encoding="utf-8")
                return LigandRecord(
                    identifier=identifier, source="smiles", value=entry["smiles"],
                    path=path, status="downloaded",
                )
            if entry.get("file"):
                source = Path(entry["file"])
                if not source.is_file():
                    raise PipelineError(f"ligand file not found: {source}")
                return LigandRecord(
                    identifier=identifier, source="file", value=str(source),
                    path=source, status="downloaded",
                )
            if entry.get("pdbqt"):
                source = Path(entry["pdbqt"])
                if not source.is_file():
                    raise PipelineError(f"ligand PDBQT not found: {source}")
                return LigandRecord(
                    identifier=identifier, source="pdbqt", value=str(source),
                    path=source, status="downloaded",
                )
            if entry.get("pubchem"):
                record = downloader.fetch_pubchem(entry["pubchem"], raw_dir)
                record.identifier = identifier  # honour the entry-level id
                return record
            if entry.get("zinc"):
                record = downloader.fetch_zinc(entry["zinc"], raw_dir)
                record.identifier = identifier
                return record
            if entry.get("pdb_ligand"):
                path = PDBDownloader().fetch_ligand(entry["pdb_ligand"], raw_dir)
                return LigandRecord(
                    identifier=identifier, source="pdb_ligand",
                    value=entry["pdb_ligand"], path=path, status="downloaded",
                )
            raise PipelineError(f"unrecognised ligand entry: {entry!r}")
        except DockFlowError as exc:
            self._log(f"{identifier}: download failed ({exc})")
            return LigandRecord(
                identifier=identifier, status="error", error=str(exc)
            )

    def _define_gridbox(self, target_record, ligand_records) -> GridBox:
        self._check_cancel()
        self._step("gridbox", "running")
        start = time.perf_counter()
        cfg = self.config.gridbox or {}
        source = cfg.get("source", "auto")
        padding = float(cfg.get("padding", 4.0))
        box: GridBox | None = None
        if cfg.get("center") and cfg.get("size"):
            box = GridBox(
                center=tuple(cfg["center"]),
                size=tuple(cfg["size"]),
                source="explicit",
                padding=padding,
            )
            source = "explicit"
        if box is None and source in ("auto", "ligand"):
            resname = cfg.get("reference_ligand_resname")
            if not resname:
                codes = target_record.ligand_codes or []
                resname = codes[0] if codes else None
            if resname and target_record.path is not None:
                try:
                    box = box_from_pocket(
                        target_record.path, resname,
                        chain=cfg.get("chain"), padding=padding,
                    )
                    source = f"ligand:{resname}"
                except (ValueError, DockFlowError) as exc:
                    logger.debug("pocket box failed: %s", exc)
        if box is None and source in ("auto", "residues") and cfg.get("residues"):
            box = box_from_residues(
                target_record.path if target_record.path else [],
                chain=cfg.get("chain"),
                residues=cfg["residues"],
                padding=padding,
            )
            source = "residues"
        if box is None and source == "auto":
            # last resort: the whole structure
            assert target_record.path is not None
            box = box_from_structure(target_record.path, padding=padding)
            source = "structure (fallback)"
        if box is None:
            raise PipelineError(
                f"could not derive a grid box (source={source!r}); set gridbox.center/size"
            )
        # Bounds checks (audit item 12): refuse impossible boxes, warn on
        # suspicious ones, and record everything in the manifest.
        box_warnings: list[str] = []
        if source.startswith("structure"):
            box_warnings.append(
                "grid box fell back to the WHOLE STRUCTURE (no reference "
                "ligand/residues found): the search space is almost certainly "
                "too large for meaningful docking; set gridbox.source or center/size"
            )
        # Box provenance (peer item 15): how strong is the pocket assumption,
        # and how much of the receptor does the box swallow.
        strength = assumption_strength(source)
        if strength in ("weak", "medium"):
            box_warnings.append(
                f"Box derived from {source!r} (assumption strength: "
                f"{strength}); consider confirming the pocket with "
                "conservation analysis or a known active-site residue list."
            )
        receptor_coverage: float | None = None
        prepared = self._receptor_result
        if prepared is not None and prepared.pdbqt_path is not None:
            try:
                from .pdbio import parse_pdbqt

                rec_atoms = [a for a in parse_pdbqt(prepared.pdbqt_path).atoms
                             if not a.is_hydrogen]
                if rec_atoms:
                    inside = sum(1 for a in rec_atoms if box.contains(a.xyz))
                    receptor_coverage = 100.0 * inside / len(rec_atoms)
                    if receptor_coverage > 60.0:
                        box_warnings.append(
                            f"Box covers {receptor_coverage:.0f}% of the "
                            "receptor (heavy atoms inside the box) - likely "
                            "too large; consider a tighter pocket definition."
                        )
            except Exception:  # noqa: BLE001 - coverage is informational only
                logger.debug("receptor coverage not computable", exc_info=True)
        ligand_heavy = next(
            (record.num_heavy_atoms for record in ligand_records
             if record.status == "prepared" and record.num_heavy_atoms),
            None,
        )
        try:
            box_warnings += validate_box(box, ligand_heavy_atoms=ligand_heavy)
        except ValueError as exc:
            logger.error("%s", exc)
            raise PipelineError(str(exc)) from exc
        for warning in box_warnings:
            logger.warning("%s", warning)
            self._log(f"warning: {warning}")
        config_path = box_to_vina_config(
            box, self._run_dir / "gridbox.txt",
            extra={"exhaustiveness": self.config.docking.get("exhaustiveness", 8)},
        )
        self._report.gridbox = {
            **box.to_dict(),
            "vina_config": str(config_path),
            "warnings": box_warnings,
            # peer item 15: provenance of the pocket assumption
            "assumption_strength": strength,
            "receptor_coverage_pct": (round(receptor_coverage, 1)
                                      if receptor_coverage is not None else None),
        }
        self._report.paths["gridbox_config"] = str(config_path)
        self._log(f"grid box: {box}")
        self._step("gridbox", "done", source)
        self._report.timings["gridbox"] = round(time.perf_counter() - start, 2)
        self._progress(0.38, "grid box defined")
        return box

    def _dock(self, receptor_result, ligand_records, box: GridBox):
        self._check_cancel()
        self._step("docking", "running")
        start = time.perf_counter()
        cfg = self.config.docking or {}
        # Covalent docking (peer item 9): only backends that support a
        # reactive-atom constraint may run it; the chemistry assumption is
        # recorded, not assumed silently.
        covalent_cfg = cfg.get("covalent") or {}
        covalent_enabled = bool(covalent_cfg.get("enabled"))
        covalent_decision: Any = "off"
        covalent_residue: str | None = None
        requested_backend = cfg.get("backend", "auto")
        if covalent_enabled:
            residue = covalent_cfg.get("receptor_reactive_residue") or {}
            missing = [key for key in ("chain", "resseq", "atom")
                       if not residue.get(key)]
            if missing or not covalent_cfg.get("ligand_reactive_atom"):
                raise PipelineError(
                    "docking.covalent requires ligand_reactive_atom (index) "
                    "and receptor_reactive_residue {chain, resseq, atom}"
                )
            supported = ("gnina", "vina-reactive")
            if requested_backend not in supported:
                raise PipelineError(
                    f"covalent docking is not supported by backend "
                    f"{requested_backend!r}; supported: {', '.join(supported)} "
                    "(covalent docking needs a reactive-atom constraint that "
                    "plain Vina cannot express)"
                )
            covalent_residue = (
                f"{residue['chain']}:{int(residue['resseq'])}:"
                f"{residue['atom']}")
            bond_length = float(covalent_cfg.get("bond_length", 1.7))
            covalent_decision = {
                "enabled": True,
                "ligand_reactive_atom": covalent_cfg.get("ligand_reactive_atom"),
                "receptor_reactive_residue": residue,
                "bond_length_a": bond_length,
                "backend": requested_backend,
                "chemistry_assumption": (
                    "single covalent bond of "
                    f"{bond_length:.2f} A between ligand atom "
                    f"{covalent_cfg.get('ligand_reactive_atom')} (reactive "
                    "atom type) and "
                    f"{covalent_residue}; bond order and protonation of the "
                    "linked atoms are NOT modelled - the constraint is "
                    "geometric only"
                ),
            }
            self._log(
                f"covalent docking: reactive bond {covalent_residue} at "
                f"{bond_length:.2f} A (chemistry assumption recorded in the "
                "manifest)")
        # Flexible receptor (peer item 8): the rigid part + flex file pair.
        flex_path = (str(receptor_result.flex_pdbqt_path)
                     if receptor_result.flex_pdbqt_path else None)
        receptor_for_docking = (receptor_result.rigid_pdbqt_path
                                or receptor_result.pdbqt_path)
        vina_cfg = VinaConfig.from_gridbox(
            box,
            scoring=cfg.get("scoring", "vina"),
            exhaustiveness=int(cfg.get("exhaustiveness", 8)),
            num_modes=int(cfg.get("num_modes", 9)),
            refine=int(cfg.get("refine", 5)),
            seed=cfg.get("seed"),
            cpu=int(cfg.get("cpu", 0)),
            timeout=float(cfg.get("timeout", 3600)),
            cnn_scoring=cfg.get("cnn_scoring"),
            cnn=cfg.get("cnn"),
            flex_pdbqt=flex_path,
            covalent_residue=covalent_residue,
        )
        engine = VinaEngine(
            vina_cfg,
            backend=requested_backend,
            workdir=self._run_dir / "docking",
            vina_exec=self.app_config.vina_exec,
            smina_exec=self.app_config.smina_exec,
            gnina_exec=self.app_config.gnina_exec,
        )
        self._log(f"docking backend: {engine.backend}")
        if flex_path:
            self._log(f"flexible receptor docking: rigid part "
                      f"{Path(receptor_for_docking).name} + "
                      f"{Path(flex_path).name}")
        ligand_pdbqts = [
            record.pdbqt_path
            for record in ligand_records
            if record.status == "prepared" and record.pdbqt_path
        ]
        records = [r for r in ligand_records if r.status == "prepared"]

        # Parallelism sanity check (audit item 34): cpu x parallel must not
        # exceed the physical cores, or users double-count their budget.
        docking_warnings: list[str] = []
        cores = os.cpu_count() or 1
        cpu = int(cfg.get("cpu", 0)) or cores
        parallel = max(1, int(cfg.get("parallel", 1)))
        if cpu * parallel > cores:
            docking_warnings.append(
                f"docking.cpu x docking.parallel = {cpu} x {parallel} = "
                f"{cpu * parallel} threads but only {cores} cores available; "
                "oversubscription slows docking and makes runtimes "
                "unreliable (see USER_GUIDE: Choosing cpu and parallel)"
            )
            for warning in docking_warnings:
                logger.warning("%s", warning)
                self._log(f"warning: {warning}")

        # Checkpoint/resume (audit item 28): skip ligands that a previous
        # (interrupted) run already docked, unless --force.
        progress_path = self._run_dir / "progress.json"
        completed = self._load_progress(progress_path)
        if getattr(self, "_force", False):
            completed = {}
        elif completed:
            skipped_stems = [path.stem for path in ligand_pdbqts
                             if path.stem in completed]
            if skipped_stems:
                self._log(
                    f"resuming: {len(skipped_stems)} ligand(s) already completed "
                    f"({', '.join(skipped_stems)}); use --force to re-dock"
                )

        pending_pdbqts = [
            path for path in ligand_pdbqts
            if path.stem not in completed
        ]
        resumed_results = self._results_from_progress(completed, records)

        def progress(fraction: float, message: str) -> None:
            self._progress(0.4 + 0.35 * fraction, message)

        def on_result(result) -> None:
            """Persist progress after every ligand (crash-safe resume).

            Keyed by the output-file stem (``xk2`` for ``xk2_out.pdbqt``)
            so the resume filter below matches regardless of whether the
            ligand record identifier differs from the file stem.
            """
            key = (Path(result.out_path).stem.removesuffix("_out")
                   if result.out_path else str(result.ligand_name))
            completed[key] = {
                "ligand_name": str(result.ligand_name),
                "status": "ok" if result.ok else "error",
                "out": str(result.out_path) if result.out_path else None,
                "best_affinity": result.best_affinity,
                "runtime": round(result.runtime, 2),
                "backend": result.backend,
            }
            progress_path.write_text(
                json.dumps({"completed": completed}, indent=2), encoding="utf-8"
            )

        pending_records = [
            r for r in records
            if r.pdbqt_path is not None and Path(r.pdbqt_path).stem not in completed
        ]
        results = engine.dock_batch(
            receptor_for_docking,
            pending_pdbqts,
            out_dir=self._run_dir / "docking",
            ligand_records=pending_records,
            parallel=parallel,
            progress=progress,
            stop_event=self.stop_event,
            on_result=on_result,
        )
        results = resumed_results + results

        # Consensus scoring (peer item 16): re-dock every ligand with each
        # additional scoring function / backend and merge the ranks
        # (rank-by-rank consensus; every component is reported).
        consensus_cfg = cfg.get("consensus")
        consensus_decision: Any = "off"
        if consensus_cfg:
            if not isinstance(consensus_cfg, (list, tuple)) or not consensus_cfg:
                raise PipelineError(
                    "docking.consensus must be a non-empty list of scoring "
                    "functions, e.g. [vina, vinardo, gnina_affinity]")
            allowed = {"vina", "vinardo", "ad4", "gnina_affinity"}
            unknown = [name for name in consensus_cfg if name not in allowed]
            if unknown:
                raise PipelineError(
                    f"unknown consensus scoring {unknown}; allowed: "
                    f"{sorted(allowed)}")
            primary = cfg.get("scoring", "vina")
            extras = [name for name in consensus_cfg if name != primary]
            consensus_decision = {
                "scoring_functions": list(consensus_cfg),
                "method": "rank-by-rank (mean of per-scoring ligand ranks)",
                "note": "consensus ranks are heuristic screening statistics",
            }
            for name in extras:
                if name == "gnina_affinity":
                    if not which(self.app_config.gnina_exec or "gnina"):
                        raise PipelineError(
                            "consensus scoring 'gnina_affinity' requires the "
                            "gnina executable (not found)")
                    extra_backend = "gnina"
                else:
                    extra_backend = "python"
                extra_cfg = VinaConfig.from_gridbox(
                    box,
                    scoring=name,
                    exhaustiveness=int(cfg.get("exhaustiveness", 8)),
                    num_modes=int(cfg.get("num_modes", 9)),
                    refine=int(cfg.get("refine", 5)),
                    seed=cfg.get("seed"),
                    cpu=int(cfg.get("cpu", 0)),
                    timeout=float(cfg.get("timeout", 3600)),
                    flex_pdbqt=flex_path,
                    covalent_residue=covalent_residue,
                )
                extra_engine = VinaEngine(
                    extra_cfg, backend=extra_backend,
                    workdir=self._run_dir / "docking" / f"consensus_{name}",
                    vina_exec=self.app_config.vina_exec,
                    smina_exec=self.app_config.smina_exec,
                    gnina_exec=self.app_config.gnina_exec,
                )
                self._log(f"consensus scoring: re-docking with {name}")
                for result in results:
                    if not result.ok or result.ligand is None:
                        continue
                    try:
                        extra_result = extra_engine.dock(
                            receptor_for_docking, result.ligand.pdbqt_path,
                            out_path=(self._run_dir / "docking" /
                                      f"consensus_{name}" /
                                      f"{result.ligand_name}_out.pdbqt"),
                            ligand_record=result.ligand,
                        )
                        if extra_result.best_affinity is not None:
                            result.extra_scores[name] = extra_result.best_affinity
                    except DockFlowError as exc:
                        docking_warnings.append(
                            f"consensus scoring {name} failed for "
                            f"{result.ligand_name}: {exc}")

        summary_csv = write_summary_csv(results, self._run_dir / "docking" / "summary.csv")
        self._report.docking = {
            "backend": engine.backend,
            "num_ligands": len(ligand_pdbqts),
            "num_ok": sum(1 for r in results if r.ok),
            "resumed": len(resumed_results),
            "warnings": docking_warnings,
            "decisions": {
                "covalent": covalent_decision,
                "consensus": consensus_decision,
                "flexible_receptor": bool(flex_path),
            },
            "results": [r.to_dict() for r in results],
            "summary_csv": str(summary_csv),
        }
        self._report.paths["docking_summary"] = str(summary_csv)
        for record in ligand_records:
            if record.status == "prepared":
                record.status = "docked"
        if not any(r.ok for r in results):
            raise PipelineError("docking failed for every ligand; see logs")
        self._step("docking", "done", f"{sum(1 for r in results if r.ok)} ligands")
        self._report.timings["docking"] = round(time.perf_counter() - start, 2)
        return results

    def _load_progress(self, progress_path: Path) -> dict[str, Any]:
        """Load the progress.json checkpoint written by an interrupted run."""
        if not progress_path.is_file():
            return {}
        try:
            data = json.loads(progress_path.read_text(encoding="utf-8"))
            return data.get("completed", {})
        except (ValueError, OSError) as exc:
            logger.warning("ignoring corrupt progress.json (%s)", exc)
            return {}

    def _results_from_progress(self, completed: dict[str, Any], records):
        """Rebuild DockingResults from a checkpoint's output PDBQTs."""
        from .models import DockingResult
        from .pdbio import parse_pdbqt_results

        by_name = {record.identifier: record for record in records}
        results = []
        for name, info in completed.items():
            out = info.get("out")
            if not out or not Path(out).is_file():
                continue
            results.append(DockingResult(
                ligand=by_name.get(name),
                ligand_name=info.get("ligand_name", name),
                poses=parse_pdbqt_results(out),
                out_path=Path(out),
                runtime=float(info.get("runtime") or 0.0),
                backend=info.get("backend", "unknown"),
                error=None if info.get("status") == "ok" else "failed (resumed)",
            ))
        return results

    def _analyze(self, results, receptor_result, target_record=None, box=None) -> None:
        self._check_cancel()
        self._step("analysis", "running")
        start = time.perf_counter()
        cfg = self.config.analysis or {}
        top_poses = int(cfg.get("top_poses", 3))
        analysis_dir = self._run_dir / "analysis"
        # Redocking validation (audit item 18): when the grid box came from
        # the co-crystallized ligand, that ligand's crystal pose is the
        # ground truth - extract it and compute crystal_rmsd for every pose.
        reference_atoms = self._crystal_reference(target_record, box)
        payload: dict[str, Any] = {}
        rmsd_method: str | None = None
        for result in rank_results(results):
            if not result.ok or result.out_path is None:
                continue
            analyses = analyze_docking_result(
                result, receptor_result.pdbqt_path, top_poses=top_poses,
                reference_atoms=reference_atoms,
            )
            payload[result.ligand_name] = [a.to_dict() for a in analyses]
            # Interaction fingerprints (peer item 18): bit per residue x
            # contact type, CSV + RDKit bitvector exports.
            if analyses:
                from .analyzer import write_ifp_outputs

                write_ifp_outputs(result.ligand_name, analyses, analysis_dir)
            for analysis in analyses:
                if analysis.crystal_rmsd_method:
                    rmsd_method = analysis.crystal_rmsd_method
                if analysis.contacts:
                    from .analyzer import write_contacts_csv

                    csv_name = f"{result.ligand_name}_pose{analysis.pose_index}_contacts.csv"
                    write_contacts_csv(analysis.contacts, analysis_dir / csv_name)
        (analysis_dir / "interactions.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
        # manifest.analysis carries the pose payload plus the RMSD method
        # actually used (peer item 17).
        self._report.analysis = {
            "crystal_rmsd_method": rmsd_method,
            "poses": payload,
        }
        self._report.paths["analysis_json"] = str(analysis_dir / "interactions.json")
        self._report.paths["crystal_rmsd_reference"] = (
            f"{self._crystal_reference_resname} in {target_record.path.name}"
            if reference_atoms is not None else None
        )
        # summary.csv is rewritten now that crystal_rmsd values are known
        # (write_summary_csv reads them off the pose records), and the
        # docking results in the manifest are re-serialised so the pose
        # dicts carry the freshly computed crystal_rmsd values.
        write_summary_csv(results, self._run_dir / "docking" / "summary.csv")
        self._report.docking["results"] = [r.to_dict() for r in results]
        # Pose export to SDF (audit item 29) with SD fields for the scores
        # and the crystal RMSD - one record per pose, no new pipeline stage.
        if (self.config.docking or {}).get("export_sdf", True):
            self._export_poses_sdf(results)
        self._step("analysis", "done")
        self._report.timings["analysis"] = round(time.perf_counter() - start, 2)
        self._progress(0.82, "analysis complete")

    def _crystal_reference(self, target_record, box):
        """Extract the co-crystallized reference ligand for RMSD validation.

        The reference is used whenever the search space was derived from a
        co-crystallized ligand - detected from the resolved box source
        (``pocket:XK2``), from the explicit ``reference_ligand_resname``
        config, or from the target's ligand codes when gridbox.source is
        ``ligand``/``auto`` (covers resumed runs whose box was reloaded
        from gridbox.txt).
        """
        from .analyzer import extract_reference_atoms

        self._crystal_reference_resname = None
        if target_record is None or target_record.path is None:
            return None
        cfg = self.config.gridbox or {}
        source = getattr(box, "source", "") or ""
        resname = None
        if source.startswith(("pocket:", "ligand:")):
            # pocket:XK2(A451) -> resname XK2 (instance label optional)
            resname = source.split(":", 1)[1].split("(", 1)[0].strip()
        elif cfg.get("reference_ligand_resname"):
            resname = cfg["reference_ligand_resname"]
        elif cfg.get("source", "auto") in ("ligand", "auto") and \
                target_record.ligand_codes and cfg.get("source") == "ligand":
            resname = target_record.ligand_codes[0]
        if not resname:
            return None
        chain = cfg.get("chain")
        try:
            reference = extract_reference_atoms(target_record.path, resname, chain=chain)
            self._crystal_reference_resname = resname
            self._log(
                f"redocking validation: reference ligand {resname} "
                f"({len(reference)} heavy atoms) extracted from the crystal structure"
            )
            return reference
        except (ValueError, DockFlowError) as exc:
            logger.warning("crystal reference unavailable: %s", exc)
            return None

    def _export_poses_sdf(self, results) -> None:
        """Write ``<ligand>_poses.sdf`` (one record per pose) via Meeko."""
        from .sdf_export import export_poses_sdf

        for result in results:
            if not result.ok or result.out_path is None:
                continue
            try:
                export_poses_sdf(result, self._run_dir / "docking")
            except DockFlowError as exc:
                self._log(f"warning: SDF export failed for {result.ligand_name}: {exc}")

    def _visualize(self, results, receptor_result, box: GridBox) -> None:
        cfg = self.config.visualization or {}
        if not cfg.get("enabled", True):
            self._stages.mark_skipped("visualization", "disabled in config")
            self._step("visualization", "skipped", "disabled in config")
            return
        self._check_cancel()
        self._step("visualization", "running")
        start = time.perf_counter()
        viz_dir = self._run_dir / "visualization"
        from .visualizer import render_best_poses

        rendered: list[str] = []
        engine_choice = cfg.get("engine", "auto")
        for result in rank_results(results)[: max(1, int(cfg.get("top_poses", 5)))]:
            if not result.ok or result.out_path is None:
                continue
            try:
                images = render_best_poses(
                    receptor_result.pdbqt_path,
                    result.out_path,
                    box,
                    viz_dir,
                    result.ligand_name,
                    affinities=[p.affinity for p in result.poses],
                    engine=engine_choice,
                    session=bool(cfg.get("session", True)),
                    top=int(cfg.get("top_poses", 5)),
                    pymol_executable=self.app_config.pymol_exec,
                )
                rendered.extend(str(p) for p in images)
            except DockFlowError as exc:
                self._log(f"visualization failed for {result.ligand_name}: {exc}")
        # Interactive 3D viewer (peer item 20): always generated for a
        # completed run - self-contained HTML (3Dmol.js CDN), receptor +
        # poses + grid box + clickable contacts + crystal reference
        # overlay when redocking validation exists.
        interactive_path = None
        try:
            from .interactive_viewer import build_interactive_html

            # The viewer reads run_dir/manifest.json, which is normally
            # published at the END of the run; snapshot the in-memory
            # report now so the viewer sees the current state (the report
            # stage rewrites the file with the final, complete data).
            self._report.stages = self._stages.to_list()
            self._report.save(self._run_dir / "manifest.json")
            interactive_path = build_interactive_html(self._run_dir)
        except Exception as exc:  # noqa: BLE001 - viewer must never fail a run
            logger.warning("interactive viewer failed: %s", exc)
        self._report.visualization = {
            "images": rendered,
            "interactive_html": str(interactive_path) if interactive_path else None,
        }
        self._report.paths["interactive_html"] = (
            str(interactive_path) if interactive_path else None)
        self._report.paths["visualization_dir"] = str(viz_dir)
        self._step("visualization", "done" if rendered else "fallback",
                   f"{len(rendered)} images"
                   + (" + interactive.html" if interactive_path else ""))
        self._report.timings["visualization"] = round(time.perf_counter() - start, 2)
        self._progress(0.92, "visualization complete")

    def _catalytic_check(self, box) -> str | None:
        """Catalytic / pH-sensitive residues near this run's pocket (item 6).

        Heuristic geometry check on the cleaned receptor (full residue
        information) against the grid box; returns the warning text or
        ``None``.  Never fatal - the point is honesty, not a blockade.
        """
        try:
            from .pdbio import parse_pdb
            from .structure_qc import (
                catalytic_residue_warning,
                detect_catalytic_residues,
            )

            receptor_info = self._report.receptor or {}
            pdb_path = receptor_info.get("pdb")
            if not pdb_path or not Path(pdb_path).is_file():
                return None
            atoms = parse_pdb(pdb_path)
            findings = detect_catalytic_residues(
                atoms, box_center=box.center, box_size=box.size)
            return catalytic_residue_warning(findings)
        except Exception:  # noqa: BLE001 - informational only, never fatal
            logger.debug("catalytic-residue check failed", exc_info=True)
            return None

    def _write_report(self, results, receptor_result, box, target_record) -> None:
        self._step("report", "running")
        from .analyzer import pose_cluster_summary, pose_coordinates

        ranked = rank_results(results)
        environment = self._report.environment or scientific_environment()
        vina_version = environment.get("vina") or "unknown"
        scoring = (self.config.docking or {}).get("scoring", "vina")
        receptor_info = self._report.receptor or {}
        decisions = receptor_info.get("decisions", {})
        gridbox_info = self._report.gridbox or {}

        def _gridbox_source_text() -> str:
            """Human-readable box source (item 9): the exact ligand/residue.

            ``pocket:XK2(A263)`` becomes ``pocket:XK2 chain A residue 263``
            so a reader can look the pocket atom up in the structure.
            """
            import re

            source = str(gridbox_info.get("source", box.source) or "")
            match = re.match(r"^pocket:([\w]+)\(([\w])([\d\w]+)\)$", source)
            if match:
                return (f"pocket:{match.group(1)} chain {match.group(2)} "
                        f"residue {match.group(3)}")
            return source

        strength = gridbox_info.get("assumption_strength") or \
            assumption_strength(str(gridbox_info.get("source", box.source) or ""))
        lines: list[str] = [
            "# DockFlow-Automator run report",
            "",
            f"- run id: **{self._report.run_id}**",
            f"- target: **{target_record.identifier}** ({target_record.source})"
            + (f", {target_record.title}" if target_record.title else ""),
            f"- receptor: `{receptor_result.pdbqt_path}` ({receptor_info.get('atoms_out', '?')} atoms,"
            f" engine: {receptor_result.engine})",
            f"- grid box: center {tuple(round(v, 2) for v in box.center)},"
            f" size {tuple(round(v, 2) for v in box.size)} A"
            f" (source: {_gridbox_source_text()},"
            f" padding {gridbox_info.get('padding', 0.0)} A,"
            f" assumption: {strength})",
            f"- docking backend: {self._report.docking.get('backend', '?')}",
            "",
        ]

        # -- Degraded-mode banner (final pass, item 8) ----------------------
        # UNMISTAKABLE and above the results table: a run whose receptor
        # was prepared with the `none` engine (or any engine comparable to
        # nothing) is scientifically degraded, not "less accurate".
        banner = degraded_mode_banner(decisions, receptor_result.engine)
        if banner:
            lines += [banner, ""]

        # -- Engine comparability + metal sites (peer items 4, 10) ----------
        incomparable = decisions.get("incomparable_to") or []
        if incomparable:
            lines += [
                f"> Comparability: results from this run are NOT directly "
                f"comparable to receptors prepared with "
                f"{', '.join(incomparable)} (hydrogen/charge differences); "
                f"see manifest.receptor.decisions.comparable_to.",
                "",
            ]
        from .structure_qc import metal_warning
        metal_note = metal_warning(
            decisions.get("metal_coordination")
            if isinstance(decisions.get("metal_coordination"), list) else [],
            backend=(self.config.docking or {}).get("backend", "vina"),
        )
        if metal_note:
            lines += [f"> Metal site: {metal_note}", ""]

        # -- Catalytic residues / protonation honesty (final pass, item 6) --
        catalytic_note = self._catalytic_check(box)
        if catalytic_note:
            lines += [f"> Protonation: {catalytic_note}", ""]

        # -- Redocking banner (audit item 18) --------------------------------
        banner_pose, banner_rmsd = None, None
        for result in ranked:
            if result.ok:
                pose_rmsds = [p.crystal_rmsd for p in result.poses
                              if p.crystal_rmsd is not None]
                if pose_rmsds:
                    best = min(range(len(pose_rmsds)), key=lambda i: pose_rmsds[i])
                    banner_pose, banner_rmsd = result.ligand_name, pose_rmsds[best]
                    break
        if banner_rmsd is not None:
            verdict = "PASS" if banner_rmsd <= 2.0 else "FAIL"
            method_text = (self._report.analysis or {}).get("crystal_rmsd_method")
            method_note = (
                f" Method: {method_text}" if method_text else ""
            )
            lines += [
                f"> **Redocking pose recovery: best pose RMSD = {banner_rmsd:.2f} A"
                f" for {banner_pose}** - {verdict}",
                "",
                f"RMSD is symmetry-aware heavy-atom RMSD to the co-crystallized"
                f" ligand pose extracted from the target structure"
                f"{method_note}.",
                "The 2.0 A pass mark is a **heuristic success threshold** "
                "(community convention, not a validated cutoff for this "
                "target); see docs/interaction_criteria.md.",
                "",
            ]

        # -- Results table ----------------------------------------------------
        lines += [
            "## Results",
            "",
            "| ligand | best affinity (kcal/mol) | best crystal RMSD (A) | poses |"
            " docking score efficiency (heuristic proxy) | runtime (s) |",
            "|---|---|---|---|---|---|",
        ]
        for result in ranked:
            affinity = result.best_affinity
            affinity_text = f"{affinity:.2f}" if affinity is not None else "n/a"
            pose_rmsds = [p.crystal_rmsd for p in result.poses if p.crystal_rmsd is not None]
            rmsd_text = f"{min(pose_rmsds):.2f}" if pose_rmsds else "n/a"
            heavy = result.ligand.num_heavy_atoms if result.ligand else None
            efficiency = (f"{affinity / heavy:.3f}"
                          if affinity is not None and heavy else "n/a")
            lines.append(
                f"| {result.ligand_name} | {affinity_text} | {rmsd_text} "
                f"| {len(result.poses)} | {efficiency} | {result.runtime:.1f} |"
            )
        lines += [
            "",
            "_docking score efficiency (heuristic proxy) = affinity / "
            "heavy-atom count: a Vina-score-derived proxy, **not** "
            "experimental ligand efficiency._",
            "",
        ]

        # -- Pose clustering (audit item 20) ----------------------------------
        if ranked and ranked[0].ok and ranked[0].poses:
            best = ranked[0]
            coords = pose_coordinates(best)
            if coords:
                clusters = pose_cluster_summary(
                    coords, [p.affinity for p in best.poses], cutoff=2.0
                )
                if clusters:
                    lines += [
                        f"## Heuristic pose clustering ({best.ligand_name})",
                        "",
                        "Heuristic clustering: Kabsch RMSD, 2.0 A threshold; "
                        "no thermodynamic meaning (cluster membership says "
                        "nothing about binding free energy).",
                        "",
                        f"{len(best.poses)} poses in **{len(clusters)} cluster(s)**"
                        + (" (converged to a single binding mode)"
                           if len(clusters) == 1 else
                           " (multiple distinct binding modes)"),
                        "",
                        "| cluster | poses | mean dG (kcal/mol) | spread (A) |"
                        " representative pose |",
                        "|---|---|---|---|---|",
                    ]
                    for row in clusters:
                        lines.append(
                            f"| {row['cluster']} | {', '.join(map(str, row['poses']))} "
                            f"| {row['mean_affinity']:.2f} "
                            f"| {row['intra_cluster_rmsd_spread']:.2f} "
                            f"| {row['representative_pose']} |"
                        )
                    lines.append("")

        # -- Residue hotspots (qualified terminology, audit item 13) ----------
        lines += ["## Residue contact hotspots (best ligand, top pose)", ""]
        if ranked and ranked[0].ok and self._report.analysis:
            best = (self._report.analysis.get("poses") or {}).get(
                ranked[0].ligand_name) or []
            if best and best[0].get("residues"):
                lines += [
                    "Contact types are **geometric criteria** (distance + atom "
                    "type only, no angle criterion); see "
                    "docs/interaction_criteria.md.",
                    "",
                    "| chain | residue | geom. H-bond contacts (<=3.5 A) |"
                    " hydrophobic | ionic | metal | closest A |",
                    "|---|---|---|---|---|---|---|",
                ]
                for residue in best[0]["residues"][:15]:
                    lines.append(
                        f"| {residue['chain'] or '-'} | {residue['resname']}"
                        f"{residue['resseq']} | {residue['geom_hbond_contacts']} "
                        f"| {residue['hydrophobic_contacts']} "
                        f"| {residue['ionic_contacts']} "
                        f"| {residue['metal_contacts']} | {residue['closest']} |"
                    )
                lines.append("")

        # -- Preparation warnings (audit item 11) ------------------------------
        all_warnings = (
            [f"receptor: {w}" for w in receptor_info.get("warnings", [])]
            + [f"grid box: {w}" for w in gridbox_info.get("warnings", [])]
            + [f"docking: {w}" for w in self._report.docking.get("warnings", [])]
        )
        if all_warnings:
            lines += ["## Warnings", ""]
            lines += [f"- {warning}" for warning in all_warnings]
            lines.append("")

        # -- Assumptions footer (audit item 25) --------------------------------
        metals = decisions.get("metals_cofactors", {})
        protonation_text = str(decisions.get("protonation", "n/a"))
        if "toolkit template" in protonation_text:
            # final pass, item 6: never let "pH 7.4" read as "pKa-correct"
            protonation_text += (
                " - catalytic residues may be mis-protonated; for rigorous "
                "work, run PROPKA and pass the protonated PDB as input")
        lines += [
            "## Assumptions this run made",
            "",
            f"- receptor treated as rigid: "
            f"{'yes' if decisions.get('rigid_receptor', True) else 'no'}",
            f"- protonation: {protonation_text}",
            f"- charge model: {decisions.get('charge_model', 'gasteiger')}"
            " (Gasteiger partial charges)",
            f"- waters: {decisions.get('waters', 'removed')}",
            f"- metals/cofactors: removed {metals.get('removed', 'none')} /"
            f" kept {metals.get('kept', 'none')}",
            f"- grid box: derived from {gridbox_info.get('source', box.source)},"
            f" padding {gridbox_info.get('padding', 0.0)} A"
            f" (assumption: {strength})",
            f"- scoring: {scoring} (AutoDock Vina {vina_version}, no rescoring)",
            "- ligand tautomers/protonation: see manifest ligands[*].prep"
            " (input state by default; no enumeration unless configured)",
            "",
            "See README section *Scientific limitations and assumptions* for "
            "the general limitations behind each of these decisions, and "
            "docs/preparation_assumptions.md for the full preparation "
            "assumption table (what DockFlow does / does NOT do / what the "
            "researcher should do for publication-grade work).",
            "",
            "## Files",
            "",
            f"- manifest: `{self._run_dir / 'manifest.json'}`",
            f"- environment fingerprint: `{self._run_dir / 'environment.json'}`",
            f"- docking summary: `{self._report.paths.get('docking_summary', '-')}`",
            f"- interactions: `{self._report.paths.get('analysis_json', '-')}`",
            f"- visualization: `{self._run_dir / 'visualization'}`",
            "",
            f"_Generated by DockFlow-Automator v{__version__}_",
        ]
        report_path = self._run_dir / "report.md"
        report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self._report.paths["report"] = str(report_path)
        self._step("report", "done", str(report_path))
        self._progress(1.0, "run complete")

