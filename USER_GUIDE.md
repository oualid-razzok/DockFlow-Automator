# Using DockFlow-Automator

A practical, end-to-end guide to automated molecular docking: pick a
protein, pick your molecules, press a button, get scored poses, interaction
tables and publication-ready 3D renders — from the **GUI**, the **CLI** or a
**single YAML file**.

> Not installed yet? Follow **[BUILD_GUIDE.md](BUILD_GUIDE.md)** first
> (≈ 10 minutes), then come back. This guide assumes `dockflow info` runs.

---

## Table of contents

1. [The 30-second tour](#1-the-30-second-tour)
2. [Running the GUI application](#2-running-the-gui-application)
3. [Running the CLI](#3-running-the-cli)
4. [One-file automated runs (YAML)](#4-one-file-automated-runs-yaml)
5. [Understanding the outputs](#5-understanding-the-outputs)
6. [Choosing & interpreting results](#6-choosing--interpreting-results)
7. [Batch docking / virtual screening](#7-batch-docking--virtual-screening)
8. [Python API cookbook](#8-python-api-cookbook)
9. [Configuration reference](#9-configuration-reference)
10. [Where your files live](#10-where-your-files-live)
11. [Troubleshooting](#11-troubleshooting)
12. [FAQ](#12-faq)

---

## 1. The 30-second tour

DockFlow-Automator chains seven stages automatically:

```text
download ──► prepare receptor ──► prepare ligands ──► grid box
                                                            │
render ◄── analyze ◄── dock (AutoDock Vina) ◄───────────────┘
```

The fastest complete run (HIV-1 protease redocking — a validated example
that reproduces the known active site):

```bash
dockflow run --config examples/configs/hiv1_protease_example.yaml
```

That single command downloads PDB `1HVR` from RCSB, pulls the co-crystal
ligand XK2 plus aspirin and caffeine decoys, prepares receptor and ligands
(the modern, maintained equivalent of MGLTools' `prepare_receptor4.py` /
`prepare_ligand4.py`), derives the search box from the co-crystal ligand,
docks everything with Vina, analyzes contacts, renders PNGs and writes
`manifest.json` + `report.md`. Open `runs/hiv1_protease_redocking/report.md`
when it finishes.

## 2. Running the GUI application

### Launching

| platform | how |
|---|---|
| any (after `pip install`) | `dockflow-gui` or `dockflow gui` |
| Linux | double-click `run_dockflow.sh` (or run it in a terminal) |
| macOS | double-click `run_dockflow.command` (opens Terminal, starts the app) |
| Windows | double-click `run_dockflow.bat` |
| Docker (Linux/X11) | see `docker/Dockerfile` header |

The window is a **six-step wizard** with a progress bar of steps at the
top — you can always go back to an earlier step; nothing is destroyed.

### Step 1 — Target

- Type a **PDB id** (`1HVR`), a **UniProt accession** (`P03367` — resolved
  to the best experimental structure via PDBe), or **browse** for a local
  `.pdb` file.
- Metadata (title, resolution, chains, co-crystallized ligands) is fetched
  and displayed — the ligand names shown here are what you can pull in
  step 2 with one click.
- Tip: for AlphaFold models use their UniProt accession and pick the
  predicted structure when offered.

### Step 2 — Ligands

Add as many ligands as you like, from any mix of sources:

- **SMILES** pasted directly (3D coordinates are generated automatically),
- **local files** (`.sdf`, `.mol2`, `.pdb`, `.pdbqt`),
- **PubChem** lookups by name, CID or SMILES,
- **ZINC22** ids,
- **co-crystal ligand from the target** — one click, ideal for redocking.

### Step 3 — Prepare

- **Prepare receptor** and **Prepare ligands** are separate buttons; each
  writes `.pdbqt` files (the AutoDock exchange format) into `prepared/`.
- Options mirror the classic MGLTools flags: chain filtering, alternate
  locations (`best` occupancy or explicit), keep/drop waters and
  heteroatoms (with metal retention), hydrogens, Gasteiger charges,
  non-polar-hydrogen merging, preparation engine
  (auto → openbabel → rdkit → obabel CLI → dependency-free).
- Ligands get Meeko treatment: sanitisation, salt stripping, largest
  fragment, ETKDGv3 3D embedding, MMFF94/UFF minimisation, torsion tree.

### Step 4 — Grid box

- **From co-crystal ligand** — the box hugs the reference ligand + padding
  (the default, and what redocking papers do).
- **From residues…** — type active-site residues (e.g. `ASP25,ASP25'`).
- **Manual** — type center/size in Å.
- The **interactive 3D preview** (pure QPainter — works even without
  OpenGL) shows the receptor ribbon and the box wireframe; rotate with
  drag, zoom with the wheel. *Export* writes a Vina config file.

### Step 5 — Docking

- Set **exhaustiveness** (8 default, 16–32 for final runs), **poses**
  (num_modes), **seed** (fixed seed = reproducible), scoring
  (`vina`/`vinardo`/`ad4`), CPU count.
- Press **Start** — docking runs in background threads; the table fills in
  as ligands finish; **Cancel actually cancels** (immediately for the CLI
  backend, at the next checkpoint for the python backend).
- The engine is chosen automatically (python bindings → vina CLI → smina)
  and shown in the log pane.

### Step 6 — Results

- **Pose table**: affinity per pose (kcal/mol), RMSD to best mode.
- **Interaction table**: hydrogen bonds, hydrophobic contacts, ionic
  contacts, metal coordination, with residue hotspot counts — pick a pose
  to populate it.
- **Render PNG** (PyMOL ray-trace when installed, matplotlib otherwise),
  **Open in PyMOL** (loads the `.pse` session with receptor + poses + grid
  box), **Export CSV** (the summary table), **Open run folder**.

The bottom log pane shows every external command and event — the same log
that ends up in `logs/pipeline.log`.

## 3. Running the CLI

`dockflow --help` lists everything; each command has its own `--help`.
The six canonical commands (1HVR example):

```bash
# 1. download the target and the co-crystal ligand
dockflow download pdb     --id 1HVR --out runs/raw
dockflow download ligand  --pdb-ligand XK2 --out runs/raw
#    other ligand sources:
#    dockflow download ligand --pubchem aspirin --out runs/raw
#    dockflow download ligand --smiles "Cn1cnc2c1c(=O)n(C)c(=O)n2C" --name caffeine --out runs/raw
#    dockflow download uniprot --id P03367 --out runs/raw

# 2. prepare receptor and ligand (MGLTools-equivalent logic)
dockflow prep receptor    --in runs/raw/1hvr.pdb --out-dir runs/prepared
dockflow prep ligand      --in runs/raw/xk2.sdf --name xk2 --out-dir runs/prepared

# 3. search space from the co-crystal ligand (+4 A padding)
dockflow gridbox          --structure runs/raw/1hvr.pdb --resname XK2 --padding 4 \
                           --out runs/gridbox.txt

# 4. dock
dockflow dock             --receptor runs/prepared/receptor.pdbqt \
                           --ligands runs/prepared/xk2.pdbqt \
                           --config runs/gridbox.txt --exhaustiveness 16 \
                           --out-dir runs/docking

# 5. analyze (contacts, hotspots, clustering, efficiency)
dockflow analyze          --docking runs/docking --receptor runs/prepared/receptor.pdbqt

# 6. render
dockflow visualize        --receptor runs/prepared/receptor.pdbqt \
                           --poses runs/docking/xk2_out.pdbqt --out runs/render.png
```

Useful global flags: `--workdir`, `--verbose`, `--config-file` (persistent
options), `--dry-run` where supported. `dockflow gui` launches the desktop
app; `dockflow info` prints the environment report; `dockflow --version`
the package version. `scripts/dockflow_cli.py` runs the CLI straight from
source without installing (`python scripts/dockflow_cli.py …`).

Every stage is restartable: the downloader caches (`~/.dockflow/cache`), and
each stage only rewrites its own outputs.

## 4. One-file automated runs (YAML)

The `run` command drives the whole pipeline from a single declarative file
— the anticrcwu-style automation: everything scripted, nothing manual.

```yaml
run_id: my_screen
workdir: runs

target:
  pdb_id: 1HVR          # or uniprot: P03367 / alphafold: true / file: path.pdb

ligands:
  - id: xk2_redock
    pdb_ligand: XK2     # co-crystal ligand straight from RCSB
  - id: aspirin_decoy
    pubchem: aspirin    # name / CID / SMILES
  - id: caffeine_decoy
    smiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
  # - file: ligands/myset.sdf
  # - zinc: ZINC000000001234

receptor:
  chains: null          # e.g. [A] to keep a single chain
  keep_water: false
  keep_hetero: false
  keep_resnames: []     # e.g. [ZN] to retain catalytic metals
  altloc: best          # best-occupancy altLoc resolution
  add_hydrogens: true
  merge_nonpolar_h: true
  engine: auto          # openbabel | rdkit | openbabel-cli | none

gridbox:
  source: ligand        # ligand | residues | explicit
  reference_ligand_resname: XK2
  padding: 4.0
  # explicit alternative:
  # source: explicit
  # center: [-8.7, 15.5, 27.9]
  # size: [22, 22, 22]

docking:
  backend: auto         # python | cli | smina
  scoring: vina         # vina | vinardo | ad4
  exhaustiveness: 16
  num_modes: 9
  seed: 2026            # fixed seed => reproducible
  cpu: 0                # 0 = all cores
  parallel: 1           # concurrent ligands (CLI backend)

analysis:
  top_poses: 3
  contacts_cutoff: 5.0

visualization:
  enabled: true
  engine: auto          # pymol | matplotlib
  session: true         # also write .pse PyMOL sessions
  top_poses: 5
```

```bash
dockflow run --config my_screen.yaml
```

Progress is streamed to the console with per-stage timings; the run is
fully restartable and every artifact is logged in `manifest.json`
(config echo, file hashes, timings, results, engine/backend choices).

## 5. Understanding the outputs

Each run produces a self-contained directory:

```text
runs/<run_id>/
├── manifest.json         # machine-readable summary (inputs, hashes, results)
├── report.md             # human-readable report with result tables
├── raw/                  # downloaded structures (as fetched)
├── prepared/             # receptor.pdbqt, ligand .pdbqt files, clean PDB
├── gridbox.txt           # Vina config file of the search space
├── docking/              # <ligand>_out.pdbqt, .log, summary.csv
├── analysis/             # interactions.json, per-pose contact CSVs
├── visualization/        # rendered PNGs (+ .pse sessions)
└── logs/pipeline.log     # everything the app did, command by command
```

Key files:

- **`<ligand>_out.pdbqt`** — all poses, each preceded by
  `REMARK VINA RESULT: -11.2 0.000 0.000` (affinity, RMSD-lb, RMSD-ub).
  Load directly into PyMOL/ChimeraX if you want to explore by hand.
- **`summary.csv`** — one row per pose across all ligands:
  ligand, pose, affinity, RMSD, ligand efficiency.
- **`interactions.json`** — per pose: H-bonds, hydrophobic contacts, ionic
  contacts, metal coordination with residue/atom details, plus the residue
  hotspot table.
- **`.pse` session** — opens in open-source PyMOL with receptor cartoon,
  poses, and the CGO grid-box wireframe already set up.
- **`report.md`** — the run summarized for your lab notebook
  (config + results tables + timings).

## 6. Choosing & interpreting results

- **Affinity (kcal/mol)** — Vina's predicted ΔG; more negative = better.
  Rough calibration: ≤ −7 plausible, ≤ −9 strong, ≤ −11 usually only for
  tight binders or artefacts (always sanity-check geometry).
- **Redocking sanity check** — with a co-crystal ligand, the top pose
  should land within ~2 Å RMSD of the crystal pose and contact the same
  residues (the 1HVR example hits the canonical flap residues ILE47/50,
  ALA28, ILE84 — that's what "correct" looks like).
- **Ligand efficiency** — affinity divided by heavy-atom count; compare
  across molecules of different sizes (values ≳ −0.3 are decent).
- **Pose clustering** — poses are clustered by Kabsch RMSD; a single tight
  cluster is more convincing than nine scattered poses.
- **Interaction table** — check chemistry, not just scores: a −10 pose
  clashing with the protein or making zero H-bonds in a polar pocket is
  suspect.
- **Decoys matter** — aspirin/caffeine in the example act as negative
  controls; your real hit should clearly out-score them.

## 7. Batch docking / virtual screening

For libraries beyond GUI comfort (hundreds to thousands of ligands):

```bash
# parallel batches across a whole SDF library
python scripts/batch_dock.py \
    --receptor runs/prepared/receptor.pdbqt \
    --ligands library.sdf \
    --config runs/gridbox.txt \
    --out-dir runs/screen \
    --parallel 4 --exhaustiveness 8
```

- Prepare an SDF once (`dockflow prep ligand --in library.sdf --out-dir …`
  handles multi-record files), then dock in parallel batches.
- Lower `--exhaustiveness` (4–8) for the first pass; re-dock the top 1 %
  at 16–32.
- The CLI backend with `--parallel N` runs N vina processes; the python
  backend parallelizes inside Vina (`--cpu`). Either way, `summary.csv`
  accumulates everything for spreadsheet triage.
- `DOCKFLOW_CPU` caps cores globally if you share the machine.

## 8. Python API cookbook

Everything the GUI and CLI do is available as a library:

```python
from dockflow_core.downloader import PDBDownloader
from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions
from dockflow_core.gridbox import box_from_pocket
from dockflow_core.docker_engine import VinaConfig, VinaEngine

# 1. target
record = PDBDownloader().fetch_structure("1HVR", "raw")

# 2. receptor (MGLTools-equivalent options)
opts = ReceptorPrepOptions(engine="rdkit", chains=["A"], keep_water=False)
receptor = ReceptorPreparator(opts).prepare(record.path, "prepared")

# 3. search space
box = box_from_pocket(record.path, "XK2", padding=4.0)

# 4. docking
engine = VinaEngine(VinaConfig.from_gridbox(box, exhaustiveness=16, seed=42),
                    backend="auto", workdir="docking")
result = engine.dock(receptor.pdbqt_path, "prepared/xk2.pdbqt")
print(result.best_affinity, "kcal/mol")
```

```python
from dockflow_core.analyzer import analyze_docking_result
from dockflow_core.visualizer import render_complex

for analysis in analyze_docking_result(result, receptor.pdbqt_path, top_poses=3):
    print(analysis.affinity, [(r.resname, r.resseq, r.total)
                              for r in analysis.residue_rows[:5]])

render_complex(receptor.pdbqt_path, [result.out_path], box,
               "render.png", session_path="session.pse", engine="pymol")
```

Or drive the whole workflow in one object:

```python
from dockflow_core.pipeline import DockingPipeline, PipelineConfig
cfg = PipelineConfig.from_yaml("my_screen.yaml")
report = DockingPipeline(cfg).run(progress=print)
```

## 9. Configuration reference

Two configuration layers:

**Per-run options** (GUI widgets / CLI flags / YAML) — see §4 and the
annotated `examples/configs/hiv1_protease_example.yaml`.

**User-level settings** (`~/.dockflow/config.yaml`, editable from the GUI
Settings dialog): workdir, cache dir, explicit paths to the `vina` / `smina`
/ `pymol` / `obabel` executables, CPU count, parallelism, log level.

**Environment variables** (highest priority, great for CI/servers):

| variable | effect |
|---|---|
| `DOCKFLOW_HOME` | config + cache root (default `~/.dockflow`) |
| `DOCKFLOW_WORKDIR` | default run directory |
| `DOCKFLOW_VINA` / `DOCKFLOW_SMINA` | explicit engine executable path |
| `DOCKFLOW_PYMOL` / `DOCKFLOW_OBABEL` | explicit tool executable path |
| `DOCKFLOW_CPU` | cap cores |
| `QT_QPA_PLATFORM=offscreen` | run the GUI headless (screenshots, tests) |

## 10. Where your files live

| what | where |
|---|---|
| runs, reports, renders | `<workdir>/` (default `~/.dockflow/runs/`) |
| download cache (re-downloads are free) | `~/.dockflow/cache/` |
| user settings | `~/.dockflow/config.yaml` |
| logs of the last pipeline run | `<run>/logs/pipeline.log` |

Move a run directory anywhere — it is self-contained (all inputs are copied
in, `manifest.json` records hashes).

## 11. Troubleshooting

**`no Vina backend available`** — install the engine: `pip install vina`
(Linux wheels) or the vina executable (conda bioconda / official Windows
binary — BUILD_GUIDE §4). Check with `dockflow info`.

**`meeko is required for ligand preparation`** — `pip install "dockflow-automator[prep]"`.

**Vina rejects your receptor (`Unknown tag ROOT`)** — you fed a
ligand-style PDBQT as the receptor. Receptors are prepared with
`dockflow prep receptor` and never contain ROOT/BRANCH blocks.

**Everything works but renders are schematic** — the matplotlib fallback
renderer is active. Install `pymol-open-source` (conda-forge) for
ray-traced output and `.pse` sessions.

**GUI docking freezes the interface** — it shouldn't; docking runs in
QThreads. If the *CLI backend* appears stuck, check the log pane — the
vina process may be waiting for CPU (lower `--parallel` or set
`DOCKFLOW_CPU`).

**Windows: `vina` not found** — the install script places `vina.exe`
inside the conda env; run from an activated env (`conda activate dockflow`)
or double-click `run_dockflow.bat`, or set `DOCKFLOW_VINA` to the full
path of `vina.exe`.

**A ligand gets skipped** — its preparation failed (unsanitisable SMILES,
mixed salts). Check `logs/pipeline.log`; the pipeline continues with the
remaining ligands and lists failures in `report.md`.

**Score differences between machines** — different backends/scoring or
different CPU counts with the same *seed* can still differ slightly
(parallel scheduling is non-deterministic inside Vina). For exact
reproducibility set `seed` **and** `cpu: 1`.

More build-level issues (missing libraries, compilation errors) are in
[BUILD_GUIDE.md §10](BUILD_GUIDE.md#10-troubleshooting).

## 12. FAQ

**Is this MGLTools?** No — MGLTools is unmaintained Python 2. DockFlow
reimplements its `prepare_receptor4.py` / `prepare_ligand4.py`
parameterization semantics (hydrogens, Gasteiger charges, united atoms,
AD4 typing, torsion trees) on maintained toolkits: OpenBabel, RDKit and
Meeko.

**Which docking backend should I pick?** Whichever `auto` selects.
The three backends (python bindings, vina CLI, smina) produce identical
artifacts; the python backend gives in-process progress callbacks, the CLI
backend parallelizes across ligands, smina adds custom scoring functions.

**Can I use my own scoring function?** Use `docking.backend: smina`
(presently best for custom scoring; see the Smina docs) — DockFlow parses
its output identically.

**How do I dock into AlphaFold models?** Give the UniProt accession as the
target; DockFlow resolves and downloads the predicted model. Treat
low-pLDDT regions with care — grid-box around the pocket you trust.

**Does it work offline?** After structures/ligands are downloaded (cached),
everything else — preparation, docking, analysis, matplotlib rendering —
runs fully offline.

**How do I cite this?** Cite AutoDock Vina
([Trott & Olson 2010](https://doi.org/10.1002/jcc.21334)) and Meeko
([Eberhardt et al. 2021](https://doi.org/10.1021/acs.jcim.1c00196)) —
the science is theirs; DockFlow just automates it.

---

*Related documents: [README.md](README.md) · [BUILD_GUIDE.md](BUILD_GUIDE.md) ·
[CHANGELOG.md](CHANGELOG.md)*
