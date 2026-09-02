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
import threading
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
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
    box_from_pocket,
    box_from_residues,
    box_from_structure,
    box_to_vina_config,
)
from .models import LigandRecord
from .preparator import LigandPreparator, LigandPrepOptions, ReceptorPreparator, ReceptorPrepOptions
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
    target: dict[str, Any] = field(default_factory=dict)
    receptor: dict[str, Any] = field(default_factory=dict)
    ligands: list[dict[str, Any]] = field(default_factory=list)
    gridbox: dict[str, Any] = field(default_factory=dict)
    docking: dict[str, Any] = field(default_factory=dict)
    analysis: dict[str, Any] = field(default_factory=dict)
    visualization: dict[str, Any] = field(default_factory=dict)
    paths: dict[str, Any] = field(default_factory=dict)

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
    def run(self) -> PipelineReport:
        """Execute the whole workflow and return the report."""
        report = self._report
        started = time.perf_counter()
        run_id = self.config.run_id or timestamped_run_id("run")
        report.run_id = run_id
        base = Path(self.config.workdir)
        self._run_dir = base / run_id
        for sub in ("raw", "prepared", "docking", "analysis", "visualization", "logs"):
            (self._run_dir / sub).mkdir(parents=True, exist_ok=True)
        report.run_dir = self._run_dir
        logfile = self._run_dir / "logs" / "pipeline.log"
        setup_logging("INFO", logfile)
        logger.info("DockFlow-Automator %s starting run %s", __version__, run_id)
        self._log(f"run {run_id} -> {self._run_dir}")

        try:
            target_record = self._download()
            receptor_result = self._prepare_receptor(target_record)
            ligand_records = self._prepare_ligands()
            box = self._define_gridbox(target_record, ligand_records)
            results = self._dock(receptor_result, ligand_records, box)
            self._analyze(results, receptor_result)
            self._visualize(results, receptor_result, box)
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
        manifest = self._run_dir / "manifest.json"
        report.save(manifest)
        report.paths["manifest"] = str(manifest)
        logger.info("run %s finished in %.1fs (ok=%s)", run_id,
                    report.timings["total"], report.ok)
        return report

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
        self._report.receptor = {
            "pdbqt": str(result.pdbqt_path),
            "pdb": str(result.pdb_path),
            "engine": result.engine,
            "atoms_out": result.atoms_out,
            "warnings": result.warnings,
        }
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
        preparator = LigandPreparator(LigandPrepOptions())
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
                prep = preparator.prepare(record.path, prep_dir, record.identifier)
                record.pdbqt_path = prep.pdbqt_path
                record.status = "prepared"
                record.num_rotatable_bonds = prep.num_rotatable_bonds
                records.append(record)
                self._log(f"{record.identifier}: prepared ({prep.num_atoms} atoms)")
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
        self._step("prepare_ligands", "done", f"{len(prepared)} ligands")
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
                return downloader.fetch_pubchem(entry["pubchem"], raw_dir)
            if entry.get("zinc"):
                return downloader.fetch_zinc(entry["zinc"], raw_dir)
            if entry.get("pdb_ligand"):
                path = PDBDownloader().fetch_ligand(entry["pdb_ligand"], raw_dir)
                return LigandRecord(
                    identifier=entry["pdb_ligand"].lower(), source="pdb_ligand",
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
        config_path = box_to_vina_config(
            box, self._run_dir / "gridbox.txt",
            extra={"exhaustiveness": self.config.docking.get("exhaustiveness", 8)},
        )
        self._report.gridbox = {**box.to_dict(), "vina_config": str(config_path)}
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
        vina_cfg = VinaConfig.from_gridbox(
            box,
            scoring=cfg.get("scoring", "vina"),
            exhaustiveness=int(cfg.get("exhaustiveness", 8)),
            num_modes=int(cfg.get("num_modes", 9)),
            refine=int(cfg.get("refine", 5)),
            seed=cfg.get("seed"),
            cpu=int(cfg.get("cpu", 0)),
            timeout=float(cfg.get("timeout", 3600)),
        )
        engine = VinaEngine(
            vina_cfg,
            backend=cfg.get("backend", "auto"),
            workdir=self._run_dir / "docking",
            vina_exec=self.app_config.vina_exec,
            smina_exec=self.app_config.smina_exec,
        )
        self._log(f"docking backend: {engine.backend}")
        ligand_pdbqts = [
            record.pdbqt_path
            for record in ligand_records
            if record.status == "prepared" and record.pdbqt_path
        ]
        records = [r for r in ligand_records if r.status == "prepared"]

        def progress(fraction: float, message: str) -> None:
            self._progress(0.4 + 0.35 * fraction, message)

        results = engine.dock_batch(
            receptor_result.pdbqt_path,
            ligand_pdbqts,
            out_dir=self._run_dir / "docking",
            ligand_records=records,
            parallel=int(cfg.get("parallel", 1)),
            progress=progress,
            stop_event=self.stop_event,
        )
        summary_csv = write_summary_csv(results, self._run_dir / "docking" / "summary.csv")
        self._report.docking = {
            "backend": engine.backend,
            "num_ligands": len(ligand_pdbqts),
            "num_ok": sum(1 for r in results if r.ok),
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

    def _analyze(self, results, receptor_result) -> None:
        self._check_cancel()
        self._step("analysis", "running")
        start = time.perf_counter()
        cfg = self.config.analysis or {}
        top_poses = int(cfg.get("top_poses", 3))
        analysis_dir = self._run_dir / "analysis"
        payload: dict[str, Any] = {}
        for result in rank_results(results):
            if not result.ok or result.out_path is None:
                continue
            analyses = analyze_docking_result(
                result, receptor_result.pdbqt_path, top_poses=top_poses
            )
            payload[result.ligand_name] = [a.to_dict() for a in analyses]
            for analysis in analyses:
                if analysis.contacts:
                    from .analyzer import write_contacts_csv

                    csv_name = f"{result.ligand_name}_pose{analysis.pose_index}_contacts.csv"
                    write_contacts_csv(analysis.contacts, analysis_dir / csv_name)
        (analysis_dir / "interactions.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
        self._report.analysis = payload
        self._report.paths["analysis_json"] = str(analysis_dir / "interactions.json")
        self._step("analysis", "done")
        self._report.timings["analysis"] = round(time.perf_counter() - start, 2)
        self._progress(0.82, "analysis complete")

    def _visualize(self, results, receptor_result, box: GridBox) -> None:
        cfg = self.config.visualization or {}
        if not cfg.get("enabled", True):
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
        self._report.visualization = {"images": rendered}
        self._report.paths["visualization_dir"] = str(viz_dir)
        self._step("visualization", "done" if rendered else "fallback",
                   f"{len(rendered)} images")
        self._report.timings["visualization"] = round(time.perf_counter() - start, 2)
        self._progress(0.92, "visualization complete")

    def _write_report(self, results, receptor_result, box, target_record) -> None:
        self._step("report", "running")
        ranked = rank_results(results)
        lines: list[str] = [
            "# DockFlow-Automator run report",
            "",
            f"- run id: **{self._report.run_id}**",
            f"- target: **{target_record.identifier}** ({target_record.source})"
            + (f", {target_record.title}" if target_record.title else ""),
            f"- receptor: `{receptor_result.pdbqt_path}` ({receptor_result.atoms_out} atoms,"
            f" engine: {receptor_result.engine})",
            f"- grid box: center {tuple(round(v, 2) for v in box.center)},"
            f" size {tuple(round(v, 2) for v in box.size)} A",
            f"- docking backend: {self._report.docking.get('backend', '?')}",
            "",
            "## Results",
            "",
            "| ligand | best affinity (kcal/mol) | poses | runtime (s) |",
            "|---|---|---|---|",
        ]
        for result in ranked:
            affinity = result.best_affinity
            affinity_text = f"{affinity:.2f}" if affinity is not None else "n/a"
            lines.append(
                f"| {result.ligand_name} | {affinity_text} "
                f"| {len(result.poses)} | {result.runtime:.1f} |"
            )
        lines += ["", "## Residue interaction hotspots (best ligand)", ""]
        if ranked and ranked[0].ok and self._report.analysis:
            best = self._report.analysis.get(ranked[0].ligand_name) or []
            if best and best[0].get("residues"):
                lines += ["| chain | residue | hbonds | hydrophobic | ionic | metal | closest A |",
                          "|---|---|---|---|---|---|---|"]
                for residue in best[0]["residues"][:15]:
                    lines.append(
                        f"| {residue['chain'] or '-'} | {residue['resname']}"
                        f"{residue['resseq']} | {residue['hbonds']} "
                        f"| {residue['hydrophobic']} | {residue['ionic']} "
                        f"| {residue['metal']} | {residue['closest']} |"
                    )
        lines += [
            "",
            "## Files",
            "",
            f"- manifest: `{self._run_dir / 'manifest.json'}`",
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

