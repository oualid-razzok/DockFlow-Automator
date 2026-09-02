# DockFlow-Automator

**Unified, automated molecular docking: target/ligand download → preparation → grid box → docking → 3D visualization — end to end.**

[![build](https://github.com/dockflow/DockFlow-Automator/actions/workflows/build.yml/badge.svg)](https://github.com/dockflow/DockFlow-Automator/actions/workflows/build.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-261230.svg)](https://github.com/astral-sh/ruff)

DockFlow-Automator wraps the modern scientific-docking stack
([AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina),
[Meeko](https://github.com/forlilab/Meeko),
[RDKit](https://www.rdkit.org),
[OpenBabel](https://github.com/openbabel/openbabel) and
[open-source PyMOL](https://github.com/schrodinger/pymol-open-source))
into one reproducible workflow with three interchangeable front ends:

| front end | audience | entry point |
|---|---|---|
| **Desktop GUI** (PyQt6) | medicinal chemists, students | `dockflow-gui` |
| **CLI** | scripters, HPC users | `dockflow <command>` |
| **Python API / pipeline** | developers, automation | `DockingPipeline` |
| **Docker image** | reproducibility | `docker run dockflow-automator` |

> The automation style follows the spirit of
> [`omicscodeathon/anticrcwu`](https://github.com/omicscodeathon/anticrcwu)
> (event-driven, fully scripted workflows), while the *chemistry* follows the
> classic MGLTools logic — `prepare_receptor4.py` / `prepare_ligand4.py` —
> reimplemented on maintained toolkits (the legacy MGLTools/AD4 python-2
> stack is intentionally **not** used; Meeko + RDKit + OpenBabel replace it
> with the same parameterization semantics).

---

## Table of contents

1. [Features](#features)
2. [Quick start](#quick-start)
   - [Install](#install)
   - [GUI in 6 steps](#gui-in-6-steps)
   - [CLI in 6 commands](#cli-in-6-commands)
   - [One-file automated run](#one-file-automated-run)
3. [Pipeline stages explained](#pipeline-stages-explained)
4. [Architecture](#architecture)
5. [Repository layout](#repository-layout)
6. [C++ accelerator bindings](#c-accelerator-bindings)
7. [Docker](#docker)
8. [Python API examples](#python-api-examples)
9. [Configuration reference](#configuration-reference)
10. [Testing & development](#testing--development)
11. [Validated example: 1HVR redocking](#validated-example-1hvr-redocking)
12. [Licenses & third-party components](#licenses--third-party-components)
13. [Troubleshooting](#troubleshooting)
14. [Roadmap](#roadmap)

---

## Features

**Automated end-to-end workflow** — one YAML file (or one GUI session) drives
structure download, ligand acquisition, receptor/ligand preparation, search-space
definition, AutoDock Vina docking, interaction analysis, rendering and reporting.

- **Target acquisition**: RCSB PDB entries by id, UniProt accessions resolved
  through the PDBe API to their best experimental structure, or AlphaFold
  predicted models; local files are first-class citizens.
- **Ligand acquisition**: PubChem (name / CID / SMILES), ZINC22, RCSB chemical
  components (ideal 3D coordinates), local SDF/MOL2/PDB/PDBQT/SMILES files,
  or plain SMILES strings with on-the-fly ETKDG 3D embedding.
- **Receptor preparation** (`prepare_receptor4.py` logic, modern toolkits):
  chain selection, alternate-location resolution (best-occupancy or explicit),
  water/heteroatom filtering with metal retention, hydrogen addition,
  Gasteiger partial charges, non-polar hydrogen merging (united atoms),
  AutoDock atom typing (A/NA/OA/SA/HD/...), PDBQT output. Four interchangeable
  engines: OpenBabel bindings → RDKit → `obabel` CLI → dependency-free fallback.
- **Ligand preparation** (`prepare_ligand4.py` logic via Meeko + RDKit):
  sanitisation, salt stripping, largest-fragment selection, ETKDGv3 3D
  embedding, MMFF94/UFF minimisation, optional dimorphite-dl protonation
  states, torsion-tree PDBQT.
- **Grid box**: derived from the co-crystallized ligand, from active-site
  residues, or explicit; exported/imported as Vina config files; live 3D
  preview in the GUI (rotate/zoom, no OpenGL required).
- **Docking**: official Vina **python bindings**, the **vina CLI**, or
  **Smina** — selected automatically. Batch docking across ligands with
  progress callbacks and cancellation, `--score-only` / `--local-only` modes,
  vina / vinardo / ad4 scoring.
- **Analysis**: pose parsing, per-pose binding affinities and RMSD tables,
  geometric interaction detection (hydrogen bonds, hydrophobic contacts,
  ionic contacts, metal coordination), residue "hotspot" tables, pose
  clustering by Kabsch RMSD, ligand efficiency, CSV/JSON export.
- **Visualization**: headless PyMOL (ray-traced PNG + `.pse` sessions + CGO
  wire-frame grid box, in-process or subprocess) with an automatic matplotlib
  fallback renderer for minimal installs; "Open in PyMOL" from the GUI.
- **Reproducibility**: every run writes a `manifest.json` + markdown
  `report.md` + logs; Docker image pins the whole scientific stack;
  seeds are configurable.

## Quick start

### Install

```bash
# core + ligand/receptor preparation (Meeko + RDKit)
pip install "dockflow-automator[prep]"

# + Vina python bindings, GUI, fallback renderer
pip install "dockflow-automator[all]"

# everything including OpenBabel wheels and the C++ accelerator
pip install "dockflow-automator[all,obabel]" ./bindings
```

Extras: `[prep]` Meeko+RDKit · `[obabel]` OpenBabel bindings · `[engine]` Vina
python bindings · `[gui]` PyQt6 · `[viz]` matplotlib · `[test]` pytest/ruff ·
`[dev]` all of the previous + build.

PyMOL and the Vina CLI are best installed through conda-forge/bioconda
(via `bash scripts/install_tools.sh`, which creates a complete `dockflow`
environment, or see [Docker](#docker) for a pinned image).

```bash
micromamba create -n dockflow -c conda-forge -c bioconda \
    python=3.10 openbabel pymol-open-source autodock-vina
micromamba run -n dockflow pip install "dockflow-automator[prep,engine,gui,viz]"
```

### GUI in 6 steps

```bash
dockflow-gui
```

1. **Target** — type a PDB id (`1HVR`), a UniProt accession, or browse a local
   file; metadata (title, resolution, co-crystal ligands) is fetched and shown.
2. **Ligands** — add SMILES, files, PubChem/ZINC lookups, or pull the
   co-crystallized ligand straight from the target with one click.
3. **Prepare** — one button each for receptor and ligands; engines, chain
   filters and hydrogens are configurable.
4. **Grid box** — "From co-crystal ligand" / "From residues…" or manual
   center/size; the interactive 3D preview updates live.
5. **Docking** — set exhaustiveness/poses/seed and press *Start*; the table
   fills in as ligands finish; *Cancel* actually cancels.
6. **Results** — pose table with affinities, interaction table, "Render PNG",
   "Open in PyMOL", CSV export.

### CLI in 6 commands

```bash
dockflow download pdb     --id 1HVR --out runs/raw
dockflow download ligand  --pdb-ligand XK2 --out runs/raw
dockflow prep receptor    --in runs/raw/1hvr.pdb --out-dir runs/prepared
dockflow prep ligand      --in runs/raw/xk2.sdf --name xk2 --out-dir runs/prepared
dockflow gridbox          --structure runs/raw/1hvr.pdb --resname XK2 --padding 4 \
                           --out runs/gridbox.txt
dockflow dock             --receptor runs/prepared/receptor.pdbqt \
                           --ligands runs/prepared/xk2.pdbqt \
                           --config runs/gridbox.txt --exhaustiveness 16 \
                           --out-dir runs/docking
dockflow analyze          --docking runs/docking --receptor runs/prepared/receptor.pdbqt
dockflow visualize        --receptor runs/prepared/receptor.pdbqt \
                           --poses runs/docking/xk2_out.pdbqt --out runs/render.png
```

`dockflow info` prints a full environment report (installed modules, detected
Vina backends, chosen preparation engine).

### One-file automated run

```bash
dockflow run --config examples/configs/hiv1_protease_example.yaml
```

```yaml
target:
  pdb_id: 1HVR
ligands:
  - id: xk2_redock
    pdb_ligand: XK2        # co-crystal ligand from RCSB
  - id: aspirin_decoy
    pubchem: aspirin
  - id: caffeine_decoy
    smiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
gridbox:
  source: ligand           # from the co-crystallized ligand
  reference_ligand_resname: XK2
  padding: 4.0
docking:
  exhaustiveness: 16
  seed: 2026               # reproducible
visualization:
  enabled: true
```

Each run produces a self-contained directory:

```text
runs/<run_id>/
├── manifest.json        # machine-readable summary (config, timings, results)
├── report.md            # human-readable report with result tables
├── raw/                 # downloaded target + ligands
├── prepared/            # receptor.pdbqt, ligand PDBQTs, clean PDB
├── gridbox.txt          # Vina config of the search space
├── docking/             # *_out.pdbqt, logs, summary.csv
├── analysis/            # interactions.json, per-pose contacts CSVs
├── visualization/       # rendered PNGs (+ .pse sessions)
└── logs/pipeline.log
```

---

## Pipeline stages explained

| stage | module | what happens | MGLTools equivalent |
|---|---|---|---|
| download | `dockflow_core.downloader` | RCSB/PDBe/AlphaFold/PubChem/ZINC with retries + cache | — |
| prepare receptor | `dockflow_core.preparator` | filter (chains, altLoc, waters, hetero) → add H → Gasteiger → merge non-polar H → AD4 types → PDBQT | `prepare_receptor4.py` |
| prepare ligand | `dockflow_core.preparator` | RDKit sanitise/3D/minimise → Meeko torsion tree → PDBQT | `prepare_ligand4.py` |
| grid box | `dockflow_core.gridbox` | bounding box + padding from ligand/residues/coords; Vina config I/O | `prepare_gpf`/autogrid notions |
| docking | `dockflow_core.docker_engine` | Vina python / vina CLI / smina; batch, progress, cancel | `vina` runs |
| analysis | `dockflow_core.analyzer` | contacts, hotspots, clustering, Kabsch RMSD, efficiency | — |
| visualization | `dockflow_core.visualizer` | PyMOL `.pml` generation + headless execution, matplotlib fallback | — |
| orchestration | `dockflow_core.pipeline` | event-driven automation, manifests, reports | anticrcwu-style automation |

Receptor option mapping (defaults follow `prepare_receptor4.py`):

| MGLTools flag | DockFlow option | default |
|---|---|---|
| `-A hydrogens` | `ReceptorPrepOptions.add_hydrogens` | `True` |
| `-U nphs` | `merge_nonpolar_h` | `True` |
| `-U lps` | lone pairs dropped | always |
| `-U altloc` | `altloc` (`"best"`/`"A"`/`""`) | `"best"` |
| `-C` chains | `chains` | all |
| `-w` waters | `keep_water` | `False` |
| `-e` no charges | `charge_model="zero"` | `"gasteiger"` |

## Architecture

```text
                ┌─────────────────────────────────────────────────┐
                │                front ends                       │
                │  dockflow-gui (PyQt6)   dockflow CLI   YAML run │
                └───────────────┬─────────────────────────────────┘
                                │  threads / argparse / pipeline config
                ┌───────────────▼─────────────────────────────────┐
                │  dockflow_core.pipeline  (event-driven runner)  │
                └──┬─────────┬─────────┬─────────┬───────────┬────┘
                   │         │         │         │           │
             downloader preparator gridbox docker_engine  analyzer/visualizer
                   │         │         │         │           │
        requests  RDKit/Meeko numpy  vina(py/cli/smina) PyMOL/matplotlib
        RCSB/PDBe OpenBabel  config  AutoDock-Vina   open-source PyMOL
        PubChem  obabel-cli  files  C++ engine      .pse sessions
        ZINC/AF  fallback
                                │
                ┌───────────────▼─────────────────────────────────┐
                │  dockflow_bindings (optional C++ accelerator)  │
                │  PDBQT parse · grid box · contacts · Kabsch RMSD│
                └─────────────────────────────────────────────────┘
```

Every stage is usable standalone (see the module docstrings) and degrades
gracefully: missing optional dependencies produce informative warnings and
a working (if less accurate) fallback, never a crash.

## Repository layout

```text
DockFlow-Automator/
├── .github/workflows/build.yml     # CI: lint, tests (3 OS), bindings, docker, release
├── dockflow_core/                  # backend logic & workflow orchestration
│   ├── cli.py                      # `dockflow` command line interface
│   ├── downloader.py               # PDB/UniProt/AlphaFold + PubChem/ZINC/RCSB ligands
│   ├── preparator.py               # MGLTools-equivalent receptor/ligand preparation
│   ├── gridbox.py                  # search-space computation + Vina config I/O
│   ├── docker_engine.py            # AutoDock Vina execution & scoring (3 backends)
│   ├── analyzer.py                 # interactions, RMSD, clustering, efficiency
│   ├── visualizer.py               # PyMOL (headless) + matplotlib rendering
│   ├── pipeline.py                 # end-to-end automated pipeline
│   ├── pdbio.py                    # dependency-free PDB/PDBQT reader/writer
│   ├── models.py                   # data models (records, poses, contacts)
│   ├── config.py                   # user configuration (~/.dockflow)
│   └── utils.py                    # logging, subprocess, helpers
├── dockflow_gui/                   # PyQt6 desktop application
│   ├── main_window.py              # 6-step wizard, menus, docking workers
│   ├── widgets.py                  # step bar, ligand table, grid-box editor + 3D preview
│   ├── threads.py                  # QThread workers (progress, cancel, errors)
│   ├── app.py / __main__.py        # application bootstrap
│   └── resources.py                # stylesheet + programmatic icon
├── bindings/                       # pybind11 C++ accelerator (scikit-build-core)
│   ├── src/dockflow_bindings.cpp
│   ├── CMakeLists.txt + pyproject.toml
│   └── README.md
├── scripts/                        # standalone utilities
│   ├── dockflow_cli.py             # run the CLI from source, no install
│   ├── batch_dock.py               # parallel virtual-screening batches
│   └── install_tools.sh            # conda environment bootstrap
├── docker/                         # Dockerfile, env.yaml, entrypoint, compose
├── examples/                       # YAML pipeline + shell walkthrough
├── tests/                          # pytest suite (11 modules, 133 tests)
├── CMakeLists.txt                  # top-level superbuild
├── pyproject.toml                  # packaging, extras, tooling config
├── CHANGELOG.md · LICENSE · README.md
```

## C++ accelerator bindings

The optional `dockflow_bindings` module accelerates the hot loops of the
analyzer (identical results, verified by tests that cross-check the C++ and
NumPy implementations against each other):

| function | purpose |
|---|---|
| `parse_pdbqt_atoms(text)` | fast PDBQT ATOM/HETATM parsing |
| `grid_box(coords, padding)` | bounding box + padding |
| `pairwise_min_dist(a, b)` | per-atom minimum distances |
| `min_contacts(a, b, cutoff)` | all pairs within a cutoff |
| `direct_rmsd` / `kabsch_rmsd` | RMSD; Horn quaternion method with an in-house Jacobi eigensolver (no BLAS dependency) |
| `box_corners`, `ligand_efficiency` | helpers |

```bash
pip install ./bindings            # scikit-build-core wheel
cmake -B bindings/build -S bindings -DPython3_EXECUTABLE=$(which python)   # or plain CMake
cmake --build bindings/build -j && cmake --install bindings/build
```

## Docker

```bash
docker build -f docker/Dockerfile -t dockflow-automator .
docker run --rm dockflow-automator info
docker run --rm -v "$PWD:/data" dockflow-automator run --config /data/run.yaml
# GUI (X11): docker run --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
#                -v "$PWD:/data" dockflow-automator gui
```

The image (micromamba + conda-forge/bioconda) pins OpenBabel, open-source
PyMOL and the AutoDock Vina executable, then installs the Python package and
the C++ bindings — the whole stack in one reproducible artifact.  See
`docker/docker-compose.yml` for a bind-mounted setup.

## Python API examples

```python
from dockflow_core.downloader import PDBDownloader
from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions
from dockflow_core.gridbox import box_from_pocket
from dockflow_core.docker_engine import VinaConfig, VinaEngine

record = PDBDownloader().fetch_structure("1HVR", "raw")
receptor = ReceptorPreparator(ReceptorPrepOptions(engine="rdkit")).prepare(
    record.path, "prepared")
box = box_from_pocket(record.path, "XK2", padding=4.0)

engine = VinaEngine(VinaConfig.from_gridbox(box, exhaustiveness=16, seed=42),
                    backend="auto", workdir="docking")
result = engine.dock(receptor.pdbqt_path, "prepared/xk2.pdbqt")
print(result.best_affinity, "kcal/mol")
```

```python
from dockflow_core.analyzer import analyze_docking_result
for analysis in analyze_docking_result(result, receptor.pdbqt_path, top_poses=3):
    print(analysis.affinity, [(r.resname, r.resseq, r.total)
                              for r in analysis.residue_rows[:5]])
```

## Configuration reference

All pipeline options with their defaults are documented in
`examples/configs/hiv1_protease_example.yaml` and in the
`PipelineConfig` dataclass. Highlights:

- `receptor.engine`: `auto | openbabel | openbabel-cli | rdkit | none`
- `docking.backend`: `auto | python | cli | smina`; `scoring`: `vina | vinardo | ad4`
- `gridbox.source`: `auto | ligand | residues | explicit`
- `analysis`: `top_poses`, `contacts_cutoff`
- `visualization`: `enabled`, `engine` (`auto | pymol | matplotlib`), `session`
- GUI/tool paths live in `~/.dockflow/config.yaml` and can be overridden with
  `DOCKFLOW_HOME`, `DOCKFLOW_VINA`, `DOCKFLOW_SMINA`, `DOCKFLOW_PYMOL`,
  `DOCKFLOW_OBABEL`, `DOCKFLOW_CPU` environment variables.

## Testing & development

```bash
pip install -e ".[test,prep]"
pytest -m "not network and not gui"     # offline suite (133 tests)
pytest -m gui                           # PyQt6 offscreen smoke tests
pytest -m network                       # live-API tests (opt-in)
ruff check .                            # lint
```

CI (`.github/workflows/build.yml`) runs the lint job, the offline test matrix
(Linux/macOS/Windows × Python 3.10/3.12), the GUI smoke job, builds the C++
binding wheels for all three platforms, builds the Docker image, and publishes
sdist + wheels on `v*` tags.

## Validated example: 1HVR redocking

A fully real validation run (performed during development with the Vina python
bindings, exhaustiveness 8, seed 2026):

```text
[1] target    1HVR downloaded from RCSB (co-crystal ligand XK2)
[2] receptor  1890 -> 1848 atoms (RDKit engine, +1296 hydrogens)
[3] ligand    XK2 prepared with Meeko (84 atoms, 10 rotatable bonds)
[4] grid box  center (-8.7, 15.5, 27.9), volume ~7000 A^3
[5] docking   9 poses; best -11.20 kcal/mol
[6] analysis  contacts dominated by ILE47, ILE50, ALA28, ILE84 -
              the canonical HIV-1 protease flap/active-site residues
[7] render    xk2_render_01.png
```

The interaction profile matching the known active site is exactly the sanity
check you want from a redocking experiment.

## Licenses & third-party components

DockFlow-Automator itself is MIT licensed. It orchestrates these projects —
please review *their* licenses before redistribution, in particular the
copyleft ones (using OpenBabel through its Python bindings distributes under
its terms; the dependency-free and RDKit engines avoid it entirely):

| component | license | role |
|---|---|---|
| AutoDock Vina | Apache-2.0 | docking engine |
| Meeko | LGPL-2.1 | ligand PDBQT preparation |
| RDKit | BSD-3 | chemistry engine / fallback prep |
| OpenBabel | GPL-2.0 | conversions + receptor prep engine |
| open-source PyMOL | BSD-like (PyMOL) | visualization |
| PyQt6 | GPL-3 / commercial | desktop GUI |
| pybind11 | BSD-3 | C++ accelerator |
| requests / NumPy / matplotlib / PyYAML | BSD/MIT/Apache | infrastructure |

## Troubleshooting

- **"no Vina backend available"** — install the python bindings
  (`pip install vina`) or a CLI (`conda install -c bioconda autodock-vina`);
  check `dockflow info`.
- **"meeko is required for ligand preparation"** — `pip install "dockflow-automator[prep]"`.
- **Vina rejects the receptor ("Unknown tag ROOT")** — you are feeding a
  ligand-style PDBQT as a receptor; regenerate with
  `dockflow prep receptor` (receptor files never contain ROOT blocks).
- **Rendering falls back to matplotlib** — PyMOL is not installed; install
  `pymol-open-source` (conda-forge) for ray-traced output.
- **GUI shows a black/blank 3D preview** — the preview needs only QPainter;
  if it stays empty, no structure/ligand has been loaded yet.
- **Slow batch docking** — use `--parallel` with the CLI backend, or a higher
  `--cpu` with the python backend (Vina itself is multithreaded).

## Roadmap

- Flexible side chains (Meeko reactive/flex preparation)
- AD4 maps (`autogrid`) as a first-class scoring path
- Multi-receptor / ensemble docking and consensus scoring
- 2D ligand-interaction diagrams
- Plugin API for custom scoring functions (Smina custom scores)

---

*If DockFlow-Automator helps your research, please cite
[Trott & Olson (2010), J. Comput. Chem. 31, 455-461](https://doi.org/10.1002/jcc.21334)
(AutoDock Vina) and [Eberhardt et al. (2021), J. Chem. Inf. Model. 61, 3891-3898](https://doi.org/10.1021/acs.jcim.1c00196) (Meeko).*
