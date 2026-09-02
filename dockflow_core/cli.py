"""DockFlow-Automator command line interface.

Subcommands mirror every pipeline stage::

    dockflow download pdb     --id 1HVR --out runs/raw
    dockflow download uniprot --id P29978 --out runs/raw
    dockflow download ligand  --pubchem aspirin --out runs/raw
    dockflow prep receptor    --in target.pdb --out-dir prepared
    dockflow prep ligand      --smiles "CCO" --out-dir prepared
    dockflow gridbox          --ligand ref.pdbqt --padding 4 --out box.txt
    dockflow dock             --receptor prepared/receptor.pdbqt
                               --ligands prepared/*.pdbqt --config box.txt
    dockflow analyze          --docking runs/docking --receptor prepared/receptor.pdbqt
    dockflow visualize        --receptor prepared/receptor.pdbqt
                               --poses runs/docking/lig_out.pdbqt
    dockflow run              --config examples/configs/hiv1_protease_example.yaml
    dockflow info
    dockflow gui
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .utils import DockFlowError, VersionReport, get_logger, module_version, setup_logging


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dockflow",
        description="DockFlow-Automator: automated molecular docking, end to end.",
    )
    parser.add_argument("--version", action="version", version=f"dockflow {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    subparsers = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--workdir", default=None, help="DockFlow work directory")

    # ---------------------------------------------------------------- download
    download = subparsers.add_parser(
        "download", parents=[common], help="download targets and ligands"
    )
    download_sub = download.add_subparsers(dest="kind", required=True)
    pdb_p = download_sub.add_parser("pdb", help="download a PDB structure")
    pdb_p.add_argument("--id", required=True, help="4-character PDB id (e.g. 1HVR)")
    pdb_p.add_argument("--out", default=".", help="output directory")
    pdb_p.add_argument("--format", default="pdb", choices=["pdb", "cif"])
    pdb_p.add_argument("--force", action="store_true", help="re-download even if cached")
    uniprot_p = download_sub.add_parser("uniprot", help="resolve a UniProt accession")
    uniprot_p.add_argument("--id", required=True, help="UniProt accession (e.g. P29978)")
    uniprot_p.add_argument("--out", default=".", help="output directory")
    uniprot_p.add_argument(
        "--alphafold", action="store_true", help="force the AlphaFold model"
    )
    ligand_p = download_sub.add_parser("ligand", help="download small molecules")
    ligand_p.add_argument("--pubchem", help="PubChem name, CID or SMILES")
    ligand_p.add_argument("--zinc", help="ZINC22 identifier (ZINC000000000001)")
    ligand_p.add_argument("--pdb-ligand", help="RCSB 3-letter code (e.g. MK1)")
    ligand_p.add_argument("--out", default=".", help="output directory")

    # --------------------------------------------------------------------- prep
    prep = subparsers.add_parser("prep", parents=[common], help="prepare structures")
    prep_sub = prep.add_subparsers(dest="kind", required=True)
    rec_p = prep_sub.add_parser("receptor", help="PDB -> receptor.pdbqt")
    rec_p.add_argument("--in", dest="input", required=True, help="input PDB file")
    rec_p.add_argument("--out-dir", default="prepared", help="output directory")
    rec_p.add_argument("--basename", default="receptor")
    rec_p.add_argument("--chains", default=None, help="comma-separated chains (A,B)")
    rec_p.add_argument("--keep-water", action="store_true")
    rec_p.add_argument("--keep-hetero", action="store_true", help="keep all HETATM")
    rec_p.add_argument("--keep-resnames", default="", help="always-keep HETATM names (ZN,MG)")
    rec_p.add_argument("--altloc", default="best", help="best | A | B | '' (keep all)")
    rec_p.add_argument("--no-hydrogens", action="store_true", help="skip H addition")
    rec_p.add_argument("--charge-model", default="gasteiger",
                       choices=["gasteiger", "zero"])
    rec_p.add_argument("--engine", default="auto",
                       choices=["auto", "openbabel", "openbabel-cli", "rdkit", "none"])
    lig_p = prep_sub.add_parser("ligand", help="SMILES/SDF -> ligand.pdbqt (Meeko)")
    lig_p.add_argument("--smiles", help="SMILES string")
    lig_p.add_argument("--in", dest="input", help="input SDF/MOL2/PDB file")
    lig_p.add_argument("--name", default=None, help="ligand identifier")
    lig_p.add_argument("--out-dir", default="prepared", help="output directory")
    lig_p.add_argument("--library", action="store_true",
                       help="treat --in as a multi-record SDF library")
    lig_p.add_argument("--no-minimize", action="store_true")
    lig_p.add_argument("--no-embed", action="store_true", help="skip 3D embedding")
    lig_p.add_argument("--protonate", action="store_true", help="dimorphite-dl pH states")
    lig_p.add_argument("--seed", type=int, default=42)

    # ------------------------------------------------------------------ gridbox
    grid = subparsers.add_parser("gridbox", parents=[common], help="compute the search space")
    grid.add_argument("--ligand", help="reference ligand PDB/PDBQT file")
    grid.add_argument("--structure", help="structure PDB file (for --resname/--residues)")
    grid.add_argument("--resname", help="co-crystallized ligand code inside --structure")
    grid.add_argument("--chain", default=None)
    grid.add_argument("--residues", default=None, help="comma-separated residue numbers")
    grid.add_argument("--center", default=None, help="explicit center x,y,z")
    grid.add_argument("--size", default=None, help="explicit size x,y,z")
    grid.add_argument("--padding", type=float, default=4.0)
    grid.add_argument("--out", default="gridbox.txt", help="Vina config output")

    # --------------------------------------------------------------------- dock
    dock = subparsers.add_parser("dock", parents=[common], help="run AutoDock Vina")
    dock.add_argument("--receptor", required=True, help="receptor PDBQT")
    dock.add_argument("--ligands", nargs="+", required=True,
                      help="ligand PDBQT files (globs allowed)")
    dock.add_argument("--config", help="Vina config file with center/size")
    dock.add_argument("--center", default=None, help="center x,y,z (overrides config)")
    dock.add_argument("--size", default=None, help="size x,y,z")
    dock.add_argument("--out-dir", default="docking")
    dock.add_argument("--backend", default="auto", choices=["auto", "python", "cli", "smina"])
    dock.add_argument("--scoring", default="vina", choices=["vina", "vinardo", "ad4"])
    dock.add_argument("--exhaustiveness", type=int, default=8)
    dock.add_argument("--num-modes", type=int, default=9)
    dock.add_argument("--refine", type=int, default=5)
    dock.add_argument("--seed", type=int, default=None)
    dock.add_argument("--cpu", type=int, default=0)
    dock.add_argument("--parallel", type=int, default=1, help="concurrent ligands")
    dock.add_argument("--timeout", type=float, default=3600)
    dock.add_argument("--score-only", action="store_true")
    dock.add_argument("--local-only", action="store_true")

    # ------------------------------------------------------------------ analyze
    analyze = subparsers.add_parser("analyze", parents=[common], help="analyze results")
    analyze.add_argument("--docking", required=True, help="docking output directory or PDBQT")
    analyze.add_argument("--receptor", required=True, help="receptor PDBQT")
    analyze.add_argument("--top-poses", type=int, default=3)
    analyze.add_argument("--cutoff", type=float, default=5.0)
    analyze.add_argument("--out-dir", default=None, help="analysis output directory")

    # ---------------------------------------------------------------- visualize
    visualize = subparsers.add_parser(
        "visualize", parents=[common], help="render receptor + poses"
    )
    visualize.add_argument("--receptor", required=True, help="receptor PDBQT/PDB")
    visualize.add_argument("--poses", required=True, help="docked PDBQT (multi-pose ok)")
    visualize.add_argument("--box", default=None, help="Vina config with the box")
    visualize.add_argument("--out", default="render.png")
    visualize.add_argument("--session", default=None, help=".pse session output")
    visualize.add_argument("--engine", default="auto", choices=["auto", "pymol", "mpl"])
    visualize.add_argument("--top", type=int, default=5, help="number of poses to render")

    # ---------------------------------------------------------------------- run
    run = subparsers.add_parser("run", parents=[common], help="full automated pipeline")
    run.add_argument("--config", required=True, help="pipeline YAML config")
    run.add_argument("--run-id", default=None)

    # --------------------------------------------------------------------- info
    subparsers.add_parser("info", help="print an environment report")

    # ---------------------------------------------------------------------- gui
    subparsers.add_parser("gui", help="launch the desktop application")
    return parser


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _triple(text: str | None) -> tuple[float, float, float] | None:
    if text is None:
        return None
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 3:
        raise SystemExit(f"expected 'x,y,z' but got {text!r}")
    return tuple(float(p) for p in parts)  # type: ignore[return-value]


def _glob_ligands(patterns: Sequence[str]) -> list[Path]:
    import glob

    files: list[Path] = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern)) if any(c in pattern for c in "*?[") else [pattern]
        for match in matches:
            path = Path(match)
            if not path.is_file():
                raise SystemExit(f"ligand file not found: {path}")
            files.append(path)
    if not files:
        raise SystemExit("no ligand files given")
    return files


# ---------------------------------------------------------------------------
# Command implementations
# ---------------------------------------------------------------------------
def cmd_download(args: argparse.Namespace) -> int:
    from .downloader import LigandDownloader, PDBDownloader

    if args.kind == "pdb":
        downloader = PDBDownloader()
        record = downloader.fetch_structure(args.id, args.out, fmt=args.format,
                                             force=args.force)
        print(f"downloaded {record.identifier} -> {record.path}")
        if record.title:
            print(f"  title     : {record.title}")
        if record.resolution is not None:
            print(f"  resolution: {record.resolution:.2f} A")
        if record.ligand_codes:
            print(f"  ligands   : {', '.join(record.ligand_codes)}")
        return 0
    if args.kind == "uniprot":
        downloader = PDBDownloader()
        if args.alphafold:
            record = downloader.fetch_alphafold_model(args.id, args.out)
            print(f"AlphaFold model for {record.identifier} -> {record.path}")
            return 0
        pdb_ids = downloader.pdbs_for_uniprot(args.id)
        if not pdb_ids:
            print(f"no PDB entries for {args.id}; use --alphafold for the prediction")
            return 1
        print(f"UniProt {args.id} maps to {len(pdb_ids)} entries: {', '.join(pdb_ids[:10])}")
        record = downloader.fetch_structure(pdb_ids[0], args.out)
        print(f"downloaded {record.identifier} -> {record.path}")
        return 0
    # ligand
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.pubchem:
        record = LigandDownloader().fetch_pubchem(args.pubchem, out)
        print(f"PubChem ligand {args.pubchem} -> {record.path}")
        return 0
    if args.zinc:
        record = LigandDownloader().fetch_zinc(args.zinc, out)
        print(f"ZINC ligand {args.zinc} -> {record.path}")
        return 0
    if args.pdb_ligand:
        path = PDBDownloader().fetch_ligand(args.pdb_ligand, out)
        print(f"RCSB ligand {args.pdb_ligand} -> {path}")
        return 0
    raise SystemExit("download ligand needs --pubchem, --zinc or --pdb-ligand")


def cmd_prep(args: argparse.Namespace) -> int:
    if args.kind == "receptor":
        from .preparator import ReceptorPreparator, ReceptorPrepOptions

        options = ReceptorPrepOptions(
            chains=[c.strip() for c in args.chains.split(",")] if args.chains else None,
            keep_water=args.keep_water,
            keep_hetero=args.keep_hetero,
            keep_resnames=[r.strip() for r in args.keep_resnames.split(",") if r.strip()],
            altloc=args.altloc if args.altloc != "''" else "",
            add_hydrogens=not args.no_hydrogens,
            charge_model=args.charge_model,
            engine=args.engine,
        )
        result = ReceptorPreparator(options).prepare(args.input, args.out_dir,
                                                     basename=args.basename)
        print(f"receptor PDBQT: {result.pdbqt_path} ({result.atoms_out} atoms)")
        print(f"clean PDB    : {result.pdb_path}")
        print(f"engine       : {result.engine}")
        for warning in result.warnings:
            print(f"  warning: {warning}")
        return 0
    # ligand
    from .preparator import LigandPreparator, LigandPrepOptions

    if not args.smiles and not args.input:
        raise SystemExit("prep ligand needs --smiles or --in")
    options = LigandPrepOptions(
        embed_3d=not args.no_embed,
        minimize=not args.no_minimize,
        protonate=args.protonate,
        random_seed=args.seed,
    )
    preparator = LigandPreparator(options)
    source = args.smiles if args.smiles else Path(args.input)
    if args.library and args.input:
        results = preparator.prepare_library(args.input, args.out_dir)
        ok = sum(1 for r in results if r.ok)
        print(f"prepared {ok}/{len(results)} library records -> {args.out_dir}")
        for result in results:
            if result.error:
                print(f"  FAILED {result.identifier}: {result.error}")
        return 0 if ok else 1
    result = preparator.prepare(source, args.out_dir, identifier=args.name)
    print(f"ligand PDBQT: {result.pdbqt_path} ({result.num_atoms} atoms, "
          f"{result.num_rotatable_bonds} rotatable bonds)")
    if result.sdf_path:
        print(f"3D SDF     : {result.sdf_path}")
    for warning in result.warnings:
        print(f"  warning: {warning}")
    return 0


def cmd_gridbox(args: argparse.Namespace) -> int:
    from .gridbox import (
        GridBox,
        box_from_pocket,
        box_from_reference_ligand,
        box_from_residues,
        box_to_vina_config,
    )

    box: GridBox | None = None
    if args.ligand:
        box = box_from_reference_ligand(args.ligand, padding=args.padding)
    elif args.resname and args.structure:
        box = box_from_pocket(args.structure, args.resname, chain=args.chain,
                              padding=args.padding)
    elif args.residues and args.structure:
        residues = [int(r) for r in args.residues.split(",")]
        box = box_from_residues(args.structure, args.chain, residues,
                                padding=args.padding)
    elif args.center and args.size:
        box = GridBox(center=_triple(args.center), size=_triple(args.size),
                      source="explicit", padding=args.padding)
    else:
        raise SystemExit(
            "gridbox needs --ligand, or --structure with --resname/--residues, "
            "or explicit --center and --size"
        )
    box_to_vina_config(box, args.out)
    print(box)
    print(f"vina config: {args.out}")
    return 0


def cmd_dock(args: argparse.Namespace) -> int:
    from .docker_engine import VinaConfig, VinaEngine, write_summary_csv

    config = VinaConfig.from_vina_config_file(args.config) if args.config else VinaConfig()
    if args.center:
        config.center = _triple(args.center)
    if args.size:
        config.size = _triple(args.size)
    if config.center == (0.0, 0.0, 0.0) and config.size == (20.0, 20.0, 20.0) \
            and not args.config:
        raise SystemExit("no search space: pass --config, or --center and --size")
    config.scoring = args.scoring
    config.exhaustiveness = args.exhaustiveness
    config.num_modes = args.num_modes
    config.refine = args.refine
    config.seed = args.seed
    config.cpu = args.cpu
    config.timeout = args.timeout
    ligands = _glob_ligands(args.ligands)
    engine = VinaEngine(config, backend=args.backend, workdir=args.out_dir)
    print(f"backend: {engine.backend} | ligands: {len(ligands)} | scoring: {config.scoring}")
    mode = "dock"
    if args.score_only:
        mode = "score_only"
    elif args.local_only:
        mode = "local_only"
    if mode != "dock":
        # score-only / local-only work on a single ligand at a time
        results = []
        for ligand in ligands:
            results.append(engine.dock(args.receptor, ligand, mode=mode,
                                       on_log_line=None))
    else:
        results = engine.dock_batch(args.receptor, ligands, out_dir=args.out_dir,
                                    parallel=args.parallel,
                                    progress=_print_progress)
    summary = write_summary_csv(results, Path(args.out_dir) / "summary.csv")
    print(f"\nsummary: {summary}")
    for result in results:
        if result.error:
            print(f"  {result.ligand_name}: FAILED ({result.error})")
        elif result.best_affinity is not None:
            print(f"  {result.ligand_name}: best {result.best_affinity:.2f} kcal/mol "
                  f"({len(result.poses)} poses in {result.runtime:.1f}s)")
    return 0 if any(r.ok for r in results) else 1


def _print_progress(fraction: float, message: str) -> None:
    bar = "#" * int(40 * fraction)
    print(f"\r[{bar:<40}] {int(100 * fraction):3d}% {message}", end="", flush=True)
    if fraction >= 1.0:
        print()


def cmd_analyze(args: argparse.Namespace) -> int:
    import json

    from .analyzer import analyze_docking_result
    from .docker_engine import rank_results
    from .models import DockingResult
    from .pdbio import parse_pdbqt_results

    target = Path(args.docking)
    if not target.exists():
        raise SystemExit(f"path not found: {target}")
    pdbqts = (
        sorted(target.glob("*_out.pdbqt")) if target.is_dir() else [target]
    )
    if not pdbqts:
        raise SystemExit(f"no docked PDBQT files in {target}")
    out_dir = Path(args.out_dir) if args.out_dir else (
        target.parent / "analysis" if target.is_dir() else target.with_suffix("")
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[DockingResult] = []
    for pdbqt in pdbqts:
        poses = parse_pdbqt_results(pdbqt)
        results.append(DockingResult(ligand_name=pdbqt.stem.replace("_out", ""),
                                     poses=poses, out_path=pdbqt))
    payload: dict = {}
    for result in rank_results(results):
        analyses = analyze_docking_result(result, args.receptor,
                                          top_poses=args.top_poses,
                                          cutoff=args.cutoff)
        payload[result.ligand_name] = [a.to_dict() for a in analyses]
        for analysis in analyses:
            hbonds = analysis.num_hbonds
            print(f"{result.ligand_name} pose {analysis.pose_index}: "
                  f"{analysis.affinity:.2f} kcal/mol, {analysis.num_contacts} contacts "
                  f"({hbonds} hbonds)")
            if analysis.residue_rows:
                top = analysis.residue_rows[0]
                print(f"    hotspot: {top.resname}{top.resseq} chain {top.chain or '-'} "
                      f"({top.total} contacts, closest {top.closest:.1f} A)")
    (out_dir / "interactions.json").write_text(json.dumps(payload, indent=2),
                                               encoding="utf-8")
    print(f"analysis: {out_dir / 'interactions.json'}")
    return 0


def cmd_visualize(args: argparse.Namespace) -> int:
    from .gridbox import box_from_vina_config
    from .visualizer import render_best_poses

    box = box_from_vina_config(args.box) if args.box else None
    ligand_name = Path(args.poses).stem.replace("_out", "")
    images = render_best_poses(
        args.receptor, args.poses, box, Path(args.out).parent, ligand_name,
        affinities=None, engine=args.engine if args.engine != "mpl" else "matplotlib",
        session=bool(args.session), top=args.top,
    )
    for image in images:
        print(f"rendered: {image}")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from .pipeline import DockingPipeline, PipelineConfig, PipelineEvents

    config_path = Path(args.config)
    if not config_path.is_file():
        print(f"error: pipeline config not found: {config_path}", file=sys.stderr)
        return 2
    config = PipelineConfig.from_yaml(config_path)
    if args.run_id:
        config.run_id = args.run_id
    events = PipelineEvents(
        on_step=lambda step, status, detail: print(f"[{step}] {status} {detail or ''}"),
        on_log=lambda message: print(f"    {message}"),
        on_progress=lambda fraction, message: _print_progress(fraction, message),
    )
    pipeline = DockingPipeline(config, events=events)
    report = pipeline.run()
    print()
    if report.ok:
        print(f"run {report.run_id} completed -> {report.run_dir}")
        print(f"report  : {report.run_dir / 'report.md'}")
        print(f"manifest: {report.run_dir / 'manifest.json'}")
        return 0
    print(f"run {report.run_id} failed: {report.error}")
    print(f"logs: {report.run_dir / 'logs' / 'pipeline.log'}")
    return 2


def cmd_info(_args: argparse.Namespace) -> int:
    import platform

    report = VersionReport()
    report.add("dockflow-automator", __version__)
    report.add("python", platform.python_version(), sys.executable)
    report.add("platform", f"{platform.system()} {platform.release()}")
    for module in ("requests", "numpy", "PyYAML", "rdkit", "meeko", "openbabel",
                   "vina", "PyQt6", "matplotlib", "dimorphite_dl"):
        report.add(module, module_version(module.lower()))
    try:
        import dockflow_bindings  # type: ignore

        report.add("dockflow_bindings", getattr(dockflow_bindings, "__version__", "built"))
    except ImportError:
        report.add("dockflow_bindings", None, "C++ accelerator not built")
    from .docker_engine import detect_backends

    for backend in detect_backends():
        report.add(f"backend:{backend.name}", backend.version if backend.available else None,
                   backend.detail)
    from .preparator import select_engine

    engine = select_engine("auto")
    report.add("prep engine", engine.name)
    print(report.as_text())
    return 0


def cmd_gui(_args: argparse.Namespace) -> int:
    from dockflow_gui.app import main as gui_main

    return gui_main()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
_COMMANDS = {
    "download": cmd_download,
    "prep": cmd_prep,
    "gridbox": cmd_gridbox,
    "dock": cmd_dock,
    "analyze": cmd_analyze,
    "visualize": cmd_visualize,
    "run": cmd_run,
    "info": cmd_info,
    "gui": cmd_gui,
}


def main(argv: Sequence[str] | None = None) -> int:
    """Console entry point (``dockflow``)."""
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging("DEBUG" if args.verbose else "WARNING")
    workdir = getattr(args, "workdir", None)
    if workdir:
        import os

        os.environ.setdefault("DOCKFLOW_WORKDIR", str(workdir))
    handler = _COMMANDS.get(args.command)
    if handler is None:
        parser.print_help()
        return 1
    logger = get_logger("cli")
    try:
        return handler(args)
    except DockFlowError as exc:
        logger.error("%s", exc)
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())

