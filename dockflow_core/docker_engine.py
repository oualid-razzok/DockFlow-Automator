"""AutoDock Vina docking execution and scoring.

Despite the historical module name, this has nothing to do with *Docker
containers* - "docking engine" refers to the molecular docking engine
(AutoDock Vina / Smina).

Three interchangeable backends are supported:

* ``python`` - the official ``vina`` Python bindings (``pip install vina``);
* ``cli``    - the ``vina`` command line executable;
* ``smina``  - the popular Smina fork (CLI-compatible arguments).

All backends produce the same artefacts (multi-pose PDBQT with
``REMARK VINA RESULT`` records + a log file), which are parsed back through
:mod:`dockflow_core.pdbio` so downstream code never depends on the backend.
"""

from __future__ import annotations

import re
import time
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .gridbox import GridBox, read_vina_config
from .models import DockingResult, LigandRecord, PoseRecord
from .pdbio import parse_pdbqt_results
from .utils import DockFlowError, ensure_dir, get_logger, run_command, which

logger = get_logger("engine")

__all__ = [
    "DockingEngineError",
    "VinaConfig",
    "VinaEngine",
    "detect_backends",
    "parse_vina_log",
    "rank_results",
    "write_summary_csv",
]

_PROGRESS = Callable[[float, str], None]


class DockingEngineError(DockFlowError):
    """The docking engine failed."""


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass
class VinaConfig:
    """Parameters for one Vina docking run (mirrors the vina CLI)."""

    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    size: tuple[float, float, float] = (20.0, 20.0, 20.0)
    scoring: str = "vina"  # vina | vinardo | ad4 (ad4 requires precomputed maps)
    exhaustiveness: int = 8
    num_modes: int = 9
    refine: int = 5
    seed: int | None = None
    cpu: int = 0  # 0 -> all cores
    max_evals: float | None = None
    min_rmsd: float = 1.0
    energy_range: float = 3.0
    spacing: float = 0.375
    weight_terms: dict[str, float] = field(default_factory=dict)
    timeout: float = 3600.0

    @classmethod
    def from_gridbox(cls, box: GridBox, **kwargs) -> VinaConfig:
        return cls(center=tuple(box.center), size=tuple(box.size), **kwargs)  # type: ignore[arg-type]

    @classmethod
    def from_vina_config_file(cls, path: str | Path, **kwargs) -> VinaConfig:
        data = read_vina_config(path)
        get = lambda key, default, cast=float: cast(data[key]) if key in data else default  # noqa: E731
        cfg = cls(
            center=tuple(get(f"center_{a}", 0.0) for a in "xyz"),  # type: ignore[arg-type]
            size=tuple(get(f"size_{a}", 20.0) for a in "xyz"),  # type: ignore[arg-type]
            scoring=str(data.get("scoring", "vina")),
            exhaustiveness=get("exhaustiveness", 8, int),
            num_modes=get("num_modes", 9, int),
            refine=get("refine", 5, int),
            seed=get("seed", None, int) if "seed" in data else None,
            cpu=get("cpu", 0, int),
        )
        for key, value in kwargs.items():
            setattr(cfg, key, value)
        return cfg

    # -- serialization ------------------------------------------------------
    def to_dict(self) -> dict:
        return {
            "center": list(self.center),
            "size": list(self.size),
            "scoring": self.scoring,
            "exhaustiveness": self.exhaustiveness,
            "num_modes": self.num_modes,
            "refine": self.refine,
            "seed": self.seed,
            "cpu": self.cpu,
            "min_rmsd": self.min_rmsd,
            "energy_range": self.energy_range,
            "spacing": self.spacing,
        }

    def cli_args(self, receptor: Path, ligand: Path, out: Path, log: Path,
                 mode: str = "dock") -> list[str]:
        """Build the argument vector for the vina/smina CLI backend."""
        args = [
            "--receptor", str(receptor),
            "--ligand", str(ligand),
            "--out", str(out),
            "--log", str(log),
            "--center_x", f"{self.center[0]:.4f}",
            "--center_y", f"{self.center[1]:.4f}",
            "--center_z", f"{self.center[2]:.4f}",
            "--size_x", f"{self.size[0]:.4f}",
            "--size_y", f"{self.size[1]:.4f}",
            "--size_z", f"{self.size[2]:.4f}",
            "--scoring", self.scoring,
            "--num_modes", str(self.num_modes),
            "--min_rmsd", f"{self.min_rmsd:g}",
            "--energy_range", f"{self.energy_range:g}",
        ]
        if mode == "dock":
            args += ["--exhaustiveness", str(self.exhaustiveness)]
            if self.max_evals is not None:
                args += ["--max_evals", str(int(self.max_evals))]
        if mode == "score_only":
            args.append("--score_only")
        if mode == "local_only":
            args.append("--local_only")
        if self.seed is not None:
            args += ["--seed", str(self.seed)]
        if self.cpu:
            args += ["--cpu", str(self.cpu)]
        for key, value in self.weight_terms.items():
            args += [f"--weight_{key}", str(value)]
        return args


# ---------------------------------------------------------------------------
# Backend detection
# ---------------------------------------------------------------------------
@dataclass
class BackendReport:
    name: str
    available: bool
    version: str = ""
    detail: str = ""


def detect_backends(vina_exec: str | None = None, smina_exec: str | None = None) -> list[BackendReport]:
    """Probe every backend and return availability + version info."""
    reports: list[BackendReport] = []
    try:
        import vina  # noqa: F401

        version = getattr(vina, "__version__", "1.2.x")
        reports.append(BackendReport("python", True, str(version), "vina python bindings"))
    except ImportError:
        reports.append(BackendReport("python", False, "", "pip install vina"))
    for name, executable in (("cli", vina_exec or "vina"), ("smina", smina_exec or "smina")):
        path = which(executable)
        if path:
            version = _probe_cli_version(path)
            reports.append(BackendReport(name, True, version, path))
        else:
            reports.append(BackendReport(name, False, "", f"{executable} not on PATH"))
    return reports


def _probe_cli_version(executable: str) -> str:
    try:
        result = run_command([executable, "--version"], timeout=20)
        for line in result.stdout.splitlines():
            if "ina" in line or "Smina" in line:
                return line.strip()[:80]
        return result.stdout.splitlines()[0][:80] if result.stdout else "unknown"
    except DockFlowError:
        return "unknown"


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------
_LOG_TABLE_ROW = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*$"
)


def parse_vina_log(log_text: str) -> list[dict]:
    """Parse the Vina results table from a log file.

    The table looks like::

        mode |   affinity | dist from best mode
             | (kcal/mol) | rmsd l.b. | rmsd u.b.
        -----+------------+----------+----------
           1        -10.5          0.000      0.000
    """
    rows: list[dict] = []
    in_table = False
    for line in log_text.splitlines():
        if set(line.strip()) and set(line.strip()) <= set("-+| "):
            in_table = True
            continue
        if not in_table:
            if "mode" in line and "affinity" in line:
                continue
            continue
        match = _LOG_TABLE_ROW.match(line)
        if match:
            rows.append(
                {
                    "mode": int(match.group(1)),
                    "affinity": float(match.group(2)),
                    "rmsd_lb": float(match.group(3)),
                    "rmsd_ub": float(match.group(4)),
                }
            )
    return rows


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------
class VinaEngine:
    """Execute AutoDock Vina docking through the best available backend."""

    def __init__(
        self,
        config: VinaConfig,
        backend: str = "auto",
        workdir: str | Path = ".",
        vina_exec: str | None = None,
        smina_exec: str | None = None,
    ) -> None:
        self.config = config
        self.backend_name = self._select_backend(backend, vina_exec, smina_exec)
        self.workdir = ensure_dir(workdir)
        self._vina_exec = vina_exec
        self._smina_exec = smina_exec

    # -- backend selection ----------------------------------------------------
    @staticmethod
    def _select_backend(preferred: str, vina_exec: str | None, smina_exec: str | None) -> str:
        if preferred != "auto":
            if preferred == "python":
                try:
                    import vina  # noqa: F401
                    return "python"
                except ImportError as exc:
                    raise DockingEngineError(
                        "python backend requested but the vina package is not "
                        "installed (pip install vina)"
                    ) from exc
            executable = vina_exec if preferred == "cli" else smina_exec
            if which(executable or preferred):
                return preferred
            raise DockingEngineError(
                f"backend {preferred!r} requested but its executable was not found"
            )
        try:
            import vina  # noqa: F401
            return "python"
        except ImportError:
            pass
        if which(vina_exec or "vina"):
            return "cli"
        if which(smina_exec or "smina"):
            return "smina"
        raise DockingEngineError(
            "no Vina backend available: install the python bindings "
            "(pip install 'dockflow-automator[engine]') or the vina/smina "
            "executable (conda install -c bioconda autodock-vina)"
        )

    @property
    def backend(self) -> str:
        return self.backend_name

    # -- single ligand ---------------------------------------------------------
    def dock(
        self,
        receptor_pdbqt: str | Path,
        ligand_pdbqt: str | Path,
        out_path: str | Path | None = None,
        ligand_record: LigandRecord | None = None,
        mode: str = "dock",
        on_log_line: Callable[[str], None] | None = None,
    ) -> DockingResult:
        """Dock one ligand. ``mode`` is ``dock`` | ``score_only`` | ``local_only``."""
        receptor = Path(receptor_pdbqt)
        ligand = Path(ligand_pdbqt)
        if not receptor.is_file():
            raise DockingEngineError(f"receptor PDBQT not found: {receptor}")
        if not ligand.is_file():
            raise DockingEngineError(f"ligand PDBQT not found: {ligand}")
        out = Path(out_path) if out_path else (
            self.workdir / f"{ligand.stem}_out.pdbqt"
        )
        out.parent.mkdir(parents=True, exist_ok=True)
        log_path = out.with_suffix(".log")
        name = ligand_record.identifier if ligand_record else ligand.stem
        started = time.perf_counter()
        if self.backend_name == "python":
            log_text = self._dock_python(receptor, ligand, out, mode)
        else:
            log_text = self._dock_cli(receptor, ligand, out, log_path, mode, on_log_line)
        runtime = time.perf_counter() - started
        if log_text and not log_path.is_file():
            log_path.write_text(log_text, encoding="utf-8")
        poses = parse_pdbqt_results(out) if out.is_file() else []
        if not poses:
            poses = [
                PoseRecord(model=row["mode"], affinity=row["affinity"],
                           rmsd_lb=row["rmsd_lb"], rmsd_ub=row["rmsd_ub"])
                for row in parse_vina_log(log_text)
            ]
        result = DockingResult(
            ligand=ligand_record,
            ligand_name=name,
            poses=poses,
            out_path=out if out.is_file() else None,
            log_text=log_text,
            log_path=log_path if log_path.is_file() else None,
            runtime=runtime,
            backend=self.backend_name,
        )
        if not poses:
            result.error = "docking produced no poses (check the log)"
        logger.info(
            "docked %s with %s backend in %.1fs, best affinity %.2f",
            name, self.backend_name, runtime,
            result.best_affinity if result.best_affinity is not None else float("nan"),
        )
        return result

    # -- batch -----------------------------------------------------------------
    def dock_batch(
        self,
        receptor_pdbqt: str | Path,
        ligand_pdbqts: Sequence[str | Path],
        out_dir: str | Path | None = None,
        ligand_records: Sequence[LigandRecord | None] | None = None,
        parallel: int = 1,
        progress: _PROGRESS | None = None,
        stop_event=None,
    ) -> list[DockingResult]:
        """Dock a library of ligands, optionally in parallel threads.

        CLI/subprocess backends parallelise cleanly across ligands; the
        python backend is usually best run with ``parallel=1`` and a higher
        ``cpu`` setting (Vina itself is multithreaded).
        """
        out_dir = ensure_dir(out_dir or self.workdir / "docking")
        records = list(ligand_records) if ligand_records is not None else [None] * len(ligand_pdbqts)
        if len(records) != len(ligand_pdbqts):
            raise DockingEngineError("ligand_records length does not match ligand_pdbqts")
        jobs: list[tuple[Path, LigandRecord | None]] = [
            (Path(p), records[i]) for i, p in enumerate(ligand_pdbqts)
        ]
        results: list[DockingResult | None] = [None] * len(jobs)
        total = max(1, len(jobs))
        workers = max(1, min(parallel, len(jobs)))
        completed = 0

        def _run(index: int, lig_path: Path, record: LigandRecord | None) -> DockingResult:
            if stop_event is not None and stop_event.is_set():
                return DockingResult(ligand=record, ligand_name=lig_path.stem,
                                     error="cancelled", backend=self.backend_name)
            out_path = out_dir / f"{lig_path.stem}_out.pdbqt"
            try:
                return self.dock(receptor_pdbqt, lig_path, out_path, record)
            except DockingEngineError as exc:
                return DockingResult(ligand=record, ligand_name=lig_path.stem,
                                     error=str(exc), backend=self.backend_name)

        if workers == 1:
            for index, (lig_path, record) in enumerate(jobs):
                if stop_event is not None and stop_event.is_set():
                    break
                results[index] = _run(index, lig_path, record)
                completed += 1
                if progress:
                    progress(completed / total, f"docked {completed}/{total} ligands")
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {
                    pool.submit(_run, index, lig_path, record): index
                    for index, (lig_path, record) in enumerate(jobs)
                }
                for future in as_completed(futures):
                    index = futures[future]
                    results[index] = future.result()
                    completed += 1
                    if progress:
                        progress(completed / total, f"docked {completed}/{total} ligands")
                    if stop_event is not None and stop_event.is_set():
                        for pending in futures:
                            pending.cancel()
                        break
        final: list[DockingResult] = []
        cancelled = stop_event is not None and stop_event.is_set()
        for i, r in enumerate(results):
            if r is not None:
                final.append(r)
            else:
                final.append(DockingResult(
                    ligand_name=jobs[i][0].stem,
                    error="cancelled" if cancelled else "not executed",
                    backend=self.backend_name,
                ))
        return final

    # -- scoring only ------------------------------------------------------------
    def score_only(self, receptor_pdbqt: str | Path, ligand_pdbqt: str | Path) -> dict:
        """Rescore a pose without docking (Vina ``--score_only``)."""
        result = self.dock(receptor_pdbqt, ligand_pdbqt, mode="score_only")
        scores = {"inter": None, "intra": None, "torsions": None, "total": None}
        for line in result.log_text.splitlines():
            if "Estimated Free Energy of Binding" in line:
                match = re.search(r"(-?\d+\.?\d*)", line)
                if match:
                    scores["total"] = float(match.group(1))
        for line in result.log_text.splitlines():
            for key, label in (("inter", "Inter"), ("intra", "Intra"), ("torsions", "Torsional")):
                if label in line:
                    match = re.search(r"(-?\d+\.?\d*)", line)
                    if match:
                        scores[key] = float(match.group(1))
        return scores

    # -- backends ------------------------------------------------------------
    def _dock_python(self, receptor: Path, ligand: Path, out: Path, mode: str) -> str:
        """Run docking through the official Vina Python bindings."""
        try:
            from vina import Vina
        except ImportError as exc:  # pragma: no cover - guarded at selection
            raise DockingEngineError("vina python bindings unavailable") from exc
        cfg = self.config
        kwargs: dict = {"sf_name": cfg.scoring, "cpu": cfg.cpu or 0}
        if cfg.seed is not None:
            kwargs["seed"] = cfg.seed
        vina_obj = Vina(**kwargs)
        vina_obj.set_receptor(str(receptor))
        vina_obj.set_ligand_from_file(str(ligand))
        vina_obj.compute_vina_maps(center=list(cfg.center), box_size=list(cfg.size))
        log_lines: list[str] = []
        if mode == "score_only":
            energies = vina_obj.score()
            log_lines.append(f"score_only: {energies}")
            vina_obj.write_poses(str(out), overwrite=True)
        elif mode == "local_only":
            energies = vina_obj.optimize()
            log_lines.append(f"local_only: {energies}")
            vina_obj.write_poses(str(out), n_poses=1, overwrite=True)
        else:
            vina_obj.dock(exhaustiveness=cfg.exhaustiveness, n_poses=cfg.num_modes)
            vina_obj.write_poses(str(out), n_poses=cfg.num_modes, overwrite=True)
            try:
                energies = vina_obj.energies(n_poses=cfg.num_modes)
                log_lines.append(
                    "mode  affinity  intra  torsions (python backend energies)"
                )
                for i, row in enumerate(np.atleast_2d(energies)):
                    log_lines.append(
                        f"{i + 1:4d}  {float(row[0]):8.3f}  {float(row[1]):8.3f}  "
                        f"{float(row[3]) if row.size > 3 else 0.0:8.3f}"
                    )
            except Exception as exc:  # noqa: BLE001 - energies are optional
                logger.debug("energies() unavailable: %s", exc)
        log_text = "\n".join(log_lines)
        if out.is_file():
            remarks = out.read_text(encoding="utf-8", errors="replace")
            log_text = log_text + "\n" + _extract_vina_remarks(remarks)
        return log_text

    def _dock_cli(
        self,
        receptor: Path,
        ligand: Path,
        out: Path,
        log_path: Path,
        mode: str,
        on_log_line: Callable[[str], None] | None,
    ) -> str:
        executable = which(self._vina_exec or "vina")
        if self.backend_name == "smina":
            executable = which(self._smina_exec or "smina")
        if executable is None:
            raise DockingEngineError(
                f"{self.backend_name} executable disappeared (was available at init)"
            )
        args = [executable] + self.config.cli_args(receptor, ligand, out, log_path, mode)
        if self.backend_name == "smina" and mode == "score_only":
            args = [a for a in args if a != "--score_only"] + ["--score_only"]
        result = run_command(
            args,
            timeout=self.config.timeout,
            on_output=on_log_line,
        )
        log_text = ""
        if log_path.is_file():
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
        elif result.stdout:
            log_text = result.stdout
        if not result.ok and not out.is_file():
            tail = "\n".join(result.stdout.splitlines()[-12:])
            raise DockingEngineError(
                f"vina failed with exit code {result.returncode}:\n{tail}"
            )
        return log_text


def _extract_vina_remarks(pdbqt_text: str) -> str:
    """Pull REMARK VINA RESULT lines out of a poses file for the log."""
    lines = [
        line
        for line in pdbqt_text.splitlines()
        if line.strip().startswith("REMARK VINA RESULT")
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Result helpers
# ---------------------------------------------------------------------------
def rank_results(results: Sequence[DockingResult]) -> list[DockingResult]:
    """Sort results by best affinity (most negative first), failures last."""
    ok = [r for r in results if r.ok and r.best_affinity is not None]
    failed = [r for r in results if not (r.ok and r.best_affinity is not None)]
    ok.sort(key=lambda r: r.best_affinity or 0.0)
    return ok + failed


def write_summary_csv(results: Sequence[DockingResult], path: str | Path) -> Path:
    """Write a CSV summary (ligand, pose, affinity, RMSDs) for a batch run."""
    import csv

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["ligand", "pose", "affinity_kcal_mol", "rmsd_lb", "rmsd_ub", "runtime_s",
             "backend", "error"]
        )
        for result in results:
            if not result.poses:
                writer.writerow(
                    [result.ligand_name, "", "", "", "", f"{result.runtime:.1f}",
                     result.backend, result.error or ""]
                )
                continue
            for pose in result.poses:
                writer.writerow(
                    [result.ligand_name, pose.model, f"{pose.affinity:.3f}",
                     f"{pose.rmsd_lb:.3f}", f"{pose.rmsd_ub:.3f}",
                     f"{result.runtime:.1f}", result.backend, ""]
                )
    return target
