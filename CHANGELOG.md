# Changelog

All notable changes to this project are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-09-03

### Added

- **Windows support**: `scripts/install_tools.ps1` bootstrap (miniforge +
  conda-forge openbabel/PyMOL + automatic download of the official
  AutoDock Vina 1.2.7 Windows executable, renamed `vina.exe` so the CLI
  engine backend auto-detects it) and a `run_dockflow.bat` double-click
  launcher with `--cli` passthrough. The launcher is **self-healing**: it
  verifies Python and the dependencies, guards against the Windows Store
  python alias, and offers a guided setup on first run — `[1]` full conda
  stack or `[2]` quick pip install + automatic `vina.exe` download into
  the repository folder (auto-detected via `DOCKFLOW_VINA` on the next
  start).
- **macOS support**: `run_dockflow.command` double-click launcher (Finder),
  `install_tools.sh` now works on Apple Silicon and Intel (conda-forge +
  bioconda) and degrades gracefully when the vina python bindings cannot
  be installed from wheels (the bioconda CLI backend takes over). The
  launcher is **self-healing** like the Windows one: dependency probe and
  a guided setup menu on first run — `[1]` full conda stack or `[2]` quick
  pip install (with PEP 668 `--break-system-packages` fallback for
  Homebrew pythons) + automatic download of the official Vina 1.2.7 macOS
  binary for the running architecture (arm64 / x86_64).
- **Linux launcher**: `run_dockflow.sh` (same semantics as the other two,
  including the guided setup and the arch-aware Vina 1.2.7 linux_x86_64 /
  linux_aarch64 download).
- **CI/CD**: `.github/workflows/build.yml` — ruff lint; offline test matrix
  on Ubuntu/macOS/Windows × Python 3.10/3.12; offscreen GUI smoke tests on
  all three OS; cibuildwheel builds of the C++ accelerator (cp310–cp312,
  macOS x86_64 + arm64); Docker image build + smoke test; release job that
  attaches sdist, wheel and accelerator wheels to `v*` tags.
- **Documentation**: `BUILD_GUIDE.md` (per-OS compile/install manual: pip,
  conda, C++ bindings with MSVC/Xcode/gcc, vina bindings from source,
  Docker, PyInstaller standalone executables, verification checklist,
  troubleshooting) and `USER_GUIDE.md` (GUI walkthrough, CLI reference,
  YAML reference, output interpretation, virtual screening, Python API
  cookbook, FAQ). README links both and documents the platform matrix.

### Fixed

- Corrected the AutoDock Vina Windows download URL in
  `scripts/install_tools.ps1` and `run_dockflow.bat` — the release asset
  is named `vina_1.2.7_win.exe` (dots, not underscores); the underscore
  variant returns 404. Verified against the release endpoint.
- `pyproject.toml`: `gemmi` added to the `prep`/`all` extras — meeko ≥ 0.6
  imports gemmi at runtime without declaring it, which broke
  `import meeko` on Python 3.10+ (`ModuleNotFoundError: gemmi`).
- `pyproject.toml`: `openbabel-wheel` (Linux-only wheels) is now guarded
  with a `sys_platform == 'linux'` marker so `.[obabel]`/`.[all]` do not
  fail on Windows/macOS.
- Repository URLs (README badge, project urls, Dockerfile label) point to
  the actual repository.
- **CI green across all six jobs** (fixes the first failing `build` runs):
  - `docker image`: base tag corrected to `mambaorg/micromamba:1.5.10-jammy`
    (the previous `jammy-1.5.10` ordering does not exist on Docker Hub and
    the build died in 12 s with "not found"). `docker/env.yaml` no longer
    pip-installs the application before its sources are copied into the
    image (build-order bug); the entrypoint now points at
    `/opt/conda/envs/dockflow/bin` (was `/opt/conda/env/...`) and passes
    arguments to `shell`; added a `.dockerignore`.
  - `gui smoke` (all 3 OS): `QSettings` is exported by `PyQt6.QtCore`, not
    `PyQt6.QtWidgets` — the `MainWindow` import crashed on every platform
    (`ImportError`). The worker-thread smoke tests now keep the Qt event
    loop alive (`processEvents` polling) so queued cross-thread signals are
    delivered — a plain `worker.wait()` blocks the loop and the assertions
    saw empty lists.
  - `tests` (Windows): the fake vina fixture is now a `.bat` file on
    Windows and a shebang script on POSIX, so it goes through the real
    single-executable code path everywhere (a shebang script raises
    WinError 193 "not a valid Win32 application"). `run_command` converts
    any `OSError` (not just `FileNotFoundError`/`PermissionError`) into a
    reportable `ExternalToolError`, and `which()` recognises Windows
    executable suffixes explicitly.
  - CI actions bumped to current majors (`checkout@v5`, `setup-python@v6`,
    `upload/download-artifact@v7`, `cibuildwheel@v3.1.2`) — Node 20
    deprecation warnings gone.
- `pyproject.toml`: `pandas` added to the `prep`/`all` extras — meeko ≥ 0.8
  imports pandas (via `meeko.analysis`) without declaring it, which made
  `import meeko` fail silently on clean installs and skipped the ligand
  preparation tests in CI.
- `tests/test_gui_smoke.py` uses an explicit module-level skip instead of
  `pytest.importorskip(...)` with a bare `ImportError` (pytest ≥ 9.1
  deprecation); a broken system GL stack is reported in the skip reason.

## [0.1.0] - 2026-09-02

### Added

- **Automated pipeline** (`dockflow_core.pipeline`): end-to-end workflow
  target/ligand download, preparation, grid-box definition, docking,
  analysis, visualization and reporting, driven by a single YAML config
  or through the GUI.
- **Structure & ligand downloading** (`dockflow_core.downloader`): RCSB PDB
  entries, RCSB ligand ideal coordinates, UniProt-to-PDB mappings via the
  PDBe API, AlphaFold predicted models, PubChem (name / CID / SMILES) and
  ZINC22 ligand retrieval, with caching and retries.
- **Receptor preparation** (`dockflow_core.preparator`): modern
  re-implementation of the MGLTools `prepare_receptor4.py` logic
  (chain filtering, alternate-location handling, water/heteroatom removal,
  hydrogen addition, Gasteiger charges, non-polar hydrogen merging, AD4 atom
  typing, PDBQT writing) on top of OpenBabel (Python bindings or CLI) and
  RDKit, with a dependency-free fallback engine.
- **Ligand preparation**: MGLTools `prepare_ligand4.py` equivalent built on
  Meeko + RDKit (3D embedding with ETKDG, MMFF/UFF minimisation, optional
  protonation with dimorphite-dl), with an OpenBabel fallback.
- **Grid box utilities** (`dockflow_core.gridbox`): boxes from reference
  ligands, active-site residues or explicit coordinates; Vina config file
  import/export.
- **Docking engine** (`dockflow_core.docker_engine`): AutoDock Vina
  execution through the official Python bindings, the Vina CLI or Smina,
  with batch docking, scoring-only and local-only modes, progress
  callbacks and cancellation.
- **Analysis** (`dockflow_core.analyzer`): pose parsing, interaction
  detection (H-bonds, hydrophobic contacts, ionic contacts, metal
  coordination), pose clustering, Kabsch RMSD, ligand efficiency,
  CSV/JSON reports.
- **Visualization** (`dockflow_core.visualizer`): headless PyMOL rendering
  (ray-traced PNG + session files + CGO grid box), in-process or
  subprocess operation, and a matplotlib fallback renderer.
- **PyQt6 desktop application** (`dockflow_gui`): six-step wizard
  (target, ligands, preparation, grid box, docking, results), interactive
  3D grid-box preview widget, non-blocking worker threads, live logging.
- **C++ bindings** (`bindings`): pybind11 module with fast PDBQT parsing,
  grid-box computation, Kabsch RMSD and pairwise contact evaluation,
  built with scikit-build-core / CMake.
- **CLI** (`dockflow`): subcommands for every pipeline stage plus a fully
  automated `run` command; `dockflow info` prints an environment report.
- **Docker image** guaranteeing reproducibility of OpenBabel, PyMOL and Vina.
- **GitHub Actions CI** across Linux, macOS and Windows, including C++
  binding wheel builds and a Docker image build.
- **PyTest suite** covering I/O formats, preparation logic, engine
  orchestration, analysis and the end-to-end pipeline.
