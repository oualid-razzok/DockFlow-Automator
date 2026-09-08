# Changelog

## [0.2.1] - 2026-09-08

CI-repair release: fixes the failures observed when the 0.2.0 tree was
pushed to GitHub (test matrix, docker smoke, reproducibility workflow).

### Fixed
- **Docker smoke test exit 127**: the container entrypoint's `shell`
  passthrough forwarded the literal `shell` token to bash
  (`bash shell -c ...` → "No such file or directory"); it now shifts the
  token off first, so `docker run ... shell -c "command -v vina"` works.
- **`broken_openbabel` test fixture crashed on CI matrix cells without
  openbabel** (`ModuleNotFoundError` at setup): the fixture now injects a
  stub `openbabel` package with a valid `__spec__` (so
  `is_importable`/`find_spec` still sees it) when the real wheel is absent
  — macOS, Windows and Py3.10 cells exercise the same broken-dependency
  fall-through as a full install.
- **`test_kabsch_rmsd_matches_bindings_both_paths` hard-required the C++
  accelerator** (unbuilt in the plain test matrix → `ModuleNotFoundError`):
  it now `importorskip`s the module and runs for real in the dedicated
  `bindings-equivalence` CI job, which also executes the dual-path
  kabsch/clustering analyzer tests from now on.
- **`test_pose_cluster_summary_known_structure` asserted translated
  congruent shapes form 3 clusters**: `cluster_poses` measures
  superposition (Kabsch) RMSD, so rigid-body placement cannot separate
  poses — only differing *internal geometry* can.  The test data now uses
  three genuinely distinct conformations (arm-flip variants), and the test
  documents the semantics it verifies.
- **meeko was silently unusable on clean installs**: meeko 0.8.0 imports
  `scipy` (via `meeko.receptor_pdbqt`) without declaring it, exactly like
  its undeclared `gemmi`/`pandas` imports; without scipy every ligand
  preparation fails ("meeko is required for ligand preparation").  `scipy`
  is now pinned in the `prep` and `all` extras — this is also what made
  the `reproducibility` CI workflow (1HVR seed-determinism + scientific
  tests) fail before this release.

### Changed
- `bindings-equivalence` CI job additionally runs the analyzer tests that
  exercise both the C++ and pure-NumPy RMSD paths (`-k "kabsch or cluster"`).

## [0.2.0] - 2026-09-07

Audit/peer-review release: *correctness of the "auto" machinery,
reproducibility provenance and scientific transparency*.  The feature set
is otherwise frozen behind the scientific-validation track (README).

### Fixed (P0)
- **`auto` receptor-prep engine selection now actually falls through**:
  a broken optional dependency (e.g. openbabel-wheel 3.1.1.23, which
  removed `OBElementTable`) produced `AttributeError` crashes; every
  engine failure is now caught, reported as a WARNING naming the failed
  engine and the fallback, and the next engine is tried.  The OpenBabel
  engine itself works on current wheels again (local Z->symbol table,
  `OBResidue.GetChain` for the removed `GetChainID`).
- **`dockflow info` detects the OpenBabel wheel** via the actual import
  and reports the installed `openbabel-wheel` version (it printed "not
  installed" before).
- Grid boxes above Vina's 27000 A^3 hard cap are refused; large boxes
  (> 8000 A^3 with small ligands, any axis > 30 A) warn and are recorded
  in `manifest.gridbox.warnings`.
- `--stage` resumes keep the grid-box derivation and ligand heavy-atom
  counts; ligand `id`s from YAML are honoured for all ligand sources.

### Added - provenance & reproducibility (P1)
- `manifest.json`: `duration_s`, per-stage audit trail `stages[]` (wall
  time, ISO-8601 start/stop, exit status, resolved engine/backend),
  `environment` scientific-stack fingerprint; sidecar `environment.json`
  for cross-run diffing.
- **Redocking validation**: symmetry-tolerant heavy-atom Kabsch
  `crystal_rmsd` for every pose vs the co-crystallized ligand; emitted in
  `summary.csv`, `interactions.json` and a PASS/FAIL banner in
  `report.md` (threshold <= 2.0 A).  The 1HVR example recovers XK2 at
  **1.08 A**.
- Every silent preparation decision is now explicit: receptor
  `decisions` (protonation, hydrogens, charge model, waters,
  metals/cofactors removed+kept, altlocs) and per-ligand `prep` notes
  (input-hydrogen provenance, protonation source); removed species are
  counted and warned (metals/cofactors) and surfaced in `report.md`.
- `report.md`: pose-clustering summary table, Warnings section, and an
  "Assumptions this run made" footer.
- **GNINA backend** (`docking.backend: gnina`, CNN scoring via
  `cnn_scoring`/`cnn`; ADR-0001).
- **`dockflow enrich`**: ROC AUC, EF@1%/EF@5%, BEDROC (alpha=20,
  RDKit-verified) from summary.csv + actives list, optional ROC PNG.
- `dockflow run --stage {download,prep,gridbox,dock,analyze,visualize}`
  and `--dry-run` (resolved-config preview).
- **Checkpoint/resume**: `progress.json` written after every ligand;
  interrupted batch runs resume automatically, `--force` re-docks.
- Docked poses exported to `<ligand>_poses.sdf` (one record per pose,
  SD fields incl. `crystal_rmsd`; ADR-0002, `export_sdf` config key).
- `--log-format json`: structured one-event-per-line logs.
- CPU oversubscription warning (`cpu x parallel > cores`).
- `docs/interaction_criteria.md` (exact contact criteria + limitations),
  `docs/schemas/interactions.schema.json` (validated in CI),
  `docs/preparation_validation.md` (real cross-engine numbers),
  `benchmarks/redocking/` (24-complex suite + runner),
  `benchmarks/enrichment/` (DUD-E panel), `benchmarks/mgltools_comparison/`,
  `THIRD_PARTY_LICENSES.md` + `licenses.json`, ADRs, CONTRIBUTING.
- CI: `reproducibility` workflow (same-seed run twice, byte-compare
  summary.csv), C++/NumPy equivalence job with published delta artifact,
  `scientific` pytest marker, explicit schema-conformance step.

### Changed (scientific honesty)
- Interactions are now labelled **geometric contacts** everywhere
  (report headings, CSV/JSON keys `num_geom_hbond_contacts`,
  `geom_hbond_contacts`); documented distance-only criteria (no angle
  criterion).
- `ligand_efficiency` renamed to **`docking_score_efficiency`** (a
  Vina-score proxy, not experimental ligand efficiency); old name kept
  as a deprecated alias.
- aspirin/caffeine demoted to smoke-test controls (enrichment now uses
  DUD-E panels + `dockflow enrich`).
- Docker/conda are the documented reproducible installs; plain pip is
  labelled "quick smoke test".
- README restructured (intro-first), with *Scientific limitations and
  assumptions*, per-engine *Known scientific differences* and the
  software-tests vs scientific-validation split.

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
