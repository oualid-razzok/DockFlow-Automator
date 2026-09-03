# Building & Installing DockFlow-Automator

This guide covers **every supported way** to compile and install
DockFlow-Automator on **Linux, Windows and macOS**, from the quick
pip install to the full scientific stack, the C++ accelerator, the Docker
image and a standalone desktop executable.

> If you just want to *use* the application, see **[USER_GUIDE.md](USER_GUIDE.md)**.
> Time budget: quick install ≈ 2 min · full stack ≈ 10 min · source builds ≈ 15–30 min.

---

## Table of contents

1. [What "building" means here](#1-what-building-means-here)
2. [Platform support matrix](#2-platform-support-matrix)
3. [Prerequisites](#3-prerequisites)
4. [Path A — Quick install (pip only)](#4-path-a--quick-install-pip-only)
5. [Path B — Full scientific stack (conda, recommended)](#5-path-b--full-scientific-stack-conda-recommended)
6. [Path C — Docker (fully pinned)](#6-path-c--docker-fully-pinned)
7. [Compiling from source](#7-compiling-from-source)
8. [Building a standalone desktop executable](#8-building-a-standalone-desktop-executable)
9. [Verification checklist](#9-verification-checklist)
10. [Troubleshooting](#10-troubleshooting)
11. [What CI does](#11-what-ci-does)

---

## 1. What "building" means here

DockFlow-Automator is primarily a **Python 3.10+ application** — it does not
need to be compiled to run. Three optional pieces *can* be compiled:

| component | what it is | needed for | build system |
|---|---|---|---|
| `dockflow-automator` | pure Python package | everything | hatchling (`pip install .`) |
| `dockflow_bindings` | C++17 accelerator (pybind11) | ~20–50× faster analysis on big pose files — optional, pure-NumPy fallback always available | CMake ≥ 3.18 + scikit-build-core |
| `vina` python bindings | SWIG wrapper of the AutoDock Vina C++ engine | one of three *interchangeable* docking backends — optional, `vina` CLI / `smina` backends give identical results | SWIG + Boost (upstream) |

**Key design point:** every external tool is *optional and auto-detected*.
Whatever is missing is replaced by a working fallback, so the app never
fails to start — `dockflow info` always tells you what you have.

## 2. Platform support matrix

Verified against PyPI and the official release assets (2026-09):

| component | Linux x86_64 | Windows 10/11 | macOS (Intel/ARM) |
|---|---|---|---|
| Python package (core) | pip wheel | pip wheel | pip wheel |
| RDKit, Meeko, PyQt6, matplotlib | pip wheels | pip wheels | pip wheels |
| OpenBabel (python bindings) | pip wheel (`openbabel-wheel`) | conda-forge | conda-forge / Homebrew (`open-babel`) |
| Vina python bindings | pip wheels cp38–cp312 | **no wheels** — use CLI | wheels only ≤ cp3.9 — use CLI |
| Vina CLI executable | conda (bioconda) / distro packages | official `vina_1.2.7_win.exe` | conda (bioconda) |
| PyMOL open-source | conda-forge | conda-forge | conda-forge |
| C++ accelerator bindings | gcc/clang + CMake | MSVC 2019+ + CMake | Xcode CLT + CMake |
| Docker image | yes | Docker Desktop (WSL2) | Docker Desktop |
| GUI desktop app | PyQt6 (X11/Wayland) | PyQt6 (native) | PyQt6 (native) |

> **Windows note:** the Vina *python* bindings have no Windows wheels, so on
> Windows DockFlow uses the **CLI backend** with the official Vina executable
> (`scripts/install_tools.ps1` installs it automatically and renames it
> `vina.exe`). Results are identical — the backends produce the same
> PDBQT/CSV artifacts and are parsed by the same code.

## 3. Prerequisites

| OS | minimum | recommended for full stack |
|---|---|---|
| **Linux** | Python 3.10+ (distro or pyenv), pip, venv | [Miniforge](https://github.com/conda-forge/miniforge/releases) (micromamba), gcc/g++ or clang, CMake ≥ 3.18, git |
| **Windows** | Python 3.10+ from [python.org](https://www.python.org/downloads/) ("Add to PATH" checked), PowerShell 5+ | Miniforge (`Miniforge3-Windows-x86_64.exe`), [Visual Studio 2022 Build Tools](https://visualstudio.microsoft.com/downloads/) with "Desktop development with C++" (for the accelerator), git |
| **macOS** | Python 3.10+ (Xcode CLT: `xcode-select --install`), pip | Miniforge — pick `Miniforge3-MacOSX-x86_64.sh` (Intel) or `Miniforge3-MacOSX-arm64.sh` (Apple Silicon) — CMake via `brew install cmake` |
| **any** | ~2 GB free disk | ~6 GB for the full conda stack |

## 4. Path A — Quick install (pip only)

The fastest usable setup — download, prepare, dock (CLI backend) and render
(matplotlib fallback) all work. PyMOL/OpenBabel niceties come later.

### Linux / macOS

```bash
git clone https://github.com/oualid-razzok/DockFlow-Automator.git
cd DockFlow-Automator

python3 -m venv .venv && source .venv/bin/activate   # optional but clean
pip install ".[prep,gui,viz]"
```

### Windows (PowerShell)

```powershell
git clone https://github.com/oualid-razzok/DockFlow-Automator.git
cd DockFlow-Automator

python -m venv .venv ; .venv\Scripts\Activate.ps1     # optional but clean
pip install ".[prep,gui,viz]"
```

### Getting a docking engine with pip only

The Vina python bindings install from wheels **on Linux, Python 3.8–3.12**:

```bash
pip install ".[prep,gui,viz,engine]"     # Linux: adds the vina python backend
```

On Windows/macOS, get the **CLI executable** instead:

- **Windows** — download
  [`vina_1.2.7_win.exe`](https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe)
  from the AutoDock-Vina releases, rename it to `vina.exe`, put it anywhere on
  `PATH` (or point the `DOCKFLOW_VINA` environment variable at it).
  **Simplest:** just drop `vina.exe` in the repository root next to
  `run_dockflow.bat` — the launcher detects it and points `DOCKFLOW_VINA` at
  it automatically (and its quick setup downloads it for you).
- **macOS** — `conda install -c bioconda autodock-vina` (see Path B), or
  build the CLI from source (see §7.4).

Then launch:

```bash
dockflow info        # confirms what got detected
dockflow gui         # desktop application
```

## 5. Path B — Full scientific stack (conda, recommended)

This installs **everything at pinned versions**: OpenBabel bindings,
open-source PyMOL (ray-traced rendering + `.pse` sessions), the Vina
executable *and* the python bindings where wheels exist.

### 5.1 Linux / macOS — one command

```bash
bash scripts/install_tools.sh
```

The script:

1. finds micromamba/mamba/conda ([Miniforge](https://github.com/conda-forge/miniforge) recommended),
2. creates a `dockflow` environment with `openbabel`, `pymol-open-source`,
   `autodock-vina` (bioconda), `swig` and `pip`,
3. pip-installs the repository with `[prep,gui,viz]` and tries the vina
   python bindings (skipped gracefully when no wheel exists — the CLI
   backend takes over),
4. compiles the optional C++ accelerator if a compiler is present,
5. verifies everything and prints next steps.

Activate and run:

```bash
conda activate dockflow      # or: micromamba activate dockflow
dockflow info
dockflow-gui                 # or double-click run_dockflow.sh / run_dockflow.command
```

### 5.2 Windows — one command

```powershell
powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
```

The script:

1. finds Miniforge/miniconda/anaconda,
2. creates a `dockflow` environment with `openbabel` and `pymol-open-source`
   (conda-forge — bioconda has no Windows channel),
3. **downloads the official Vina 1.2.7 Windows executable** into the
   environment as `vina.exe` (auto-detected by the CLI backend),
4. pip-installs the repository with `[prep,gui,viz]` (the `engine` extra is
   omitted — no Windows wheels; the CLI backend is used instead),
5. tries to build the C++ accelerator with MSVC (continues without it if
   the compiler is missing),
6. verifies and prints next steps.

Then either activate normally:

```powershell
conda activate dockflow
dockflow info
dockflow-gui
```

…or simply **double-click `run_dockflow.bat`** in the repository root.

### 5.3 Manual conda setup (any OS, no scripts)

```bash
micromamba create -n dockflow -c conda-forge -c bioconda \
    python=3.10 openbabel pymol-open-source autodock-vina swig pip
micromamba run -n dockflow pip install ".[prep,gui,viz]"
micromamba run -n dockflow pip install ./bindings          # optional accelerator
```

(On Windows drop `autodock-vina` from the list and use the official
`vina.exe` binary as described in §4.)

## 6. Path C — Docker (fully pinned)

One image with the whole stack pinned by `docker/env.yaml` — the most
reproducible option, ideal for servers and HPC nodes with Apptainer.

```bash
git clone https://github.com/oualid-razzok/DockFlow-Automator.git
cd DockFlow-Automator

docker build -f docker/Dockerfile -t dockflow-automator .
docker run --rm dockflow-automator info
docker run --rm -v "$PWD:/data" dockflow-automator \
    run --config /data/examples/configs/hiv1_protease_example.yaml
```

Windows/macOS: use **Docker Desktop**; replace `$PWD` with `${PWD}` in
PowerShell or `$(pwd)` in zsh. GUI via X11 forwarding is possible on
Linux (`-e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix`), see
`docker/Dockerfile` header comments.

`docker/docker-compose.yml` provides a bind-mounted variant.

## 7. Compiling from source

### 7.1 The Python package (sdist / wheel)

```bash
python -m pip install build
python -m build                  # -> dist/dockflow_automator-0.1.0.tar.gz + .whl
pip install dist/dockflow_automator-0.1.0-py3-none-any.whl
```

Editable development install:

```bash
pip install -e ".[dev,prep,gui,viz]"
```

### 7.2 The C++ accelerator bindings

The bindings are built with **CMake ≥ 3.18 and a C++17 compiler**, wrapped
by scikit-build-core. A C++ compiler is the only true requirement — CMake
itself is fetched automatically by the build system when absent.

**Option 1 — pip (recommended, per platform):**

```bash
pip install ./bindings
```

| OS | toolchain needed |
|---|---|
| Linux | `sudo apt install build-essential` (or distro equivalent) |
| Windows | Visual Studio 2022 Build Tools → workload "Desktop development with C++" (from an *x64 Native Tools Command Prompt* or a normal shell after VS setup) |
| macOS | `xcode-select --install` |

**Option 2 — standalone CMake (drops the module next to the sources):**

```bash
cmake -B build -S . -DPython3_EXECUTABLE=$(which python3)   # Windows: use full python.exe path
cmake --build build -j
cmake --install build --prefix "$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
```

**Option 3 — the superbuild with OpenBabel from source (advanced):**

```bash
cmake -B build -S . -DDOCKFLOW_BUILD_OPENBABEL=ON    # heavy; pins OpenBabel 3.1.1
```

Verify after any option:

```bash
python -c "import dockflow_bindings as b; print(callable(b.kabsch_rmsd))"
```

### 7.3 The vina python bindings from source (advanced, optional)

Only needed if you insist on the *python* backend on Windows/macOS or
Python > 3.12 on Linux. Requires **SWIG ≥ 3, Boost ≥ 1.60 (headers +
python libs), CMake**:

```bash
conda install -c conda-forge swig boost-cpp boost-python-devel cmake
git clone --recursive https://github.com/ccsb-scripps/AutoDock-Vina.git
cd AutoDock-Vina
pip install .
```

The CLI backend (§4) is functionally equivalent and much easier — prefer it
unless you need the in-process API.

### 7.4 The AutoDock Vina CLI from source (any OS)

```bash
git clone https://github.com/ccsb-scripps/AutoDock-Vina.git
cd AutoDock-Vina && mkdir build && cd build
cmake .. && make && sudo make install     # Windows: cmake --build . --config Release
```

On Windows the prebuilt `vina_1.2.7_win.exe` (§4) makes this unnecessary.

## 8. Building a standalone desktop executable

You can freeze the GUI into a double-clickable app folder with PyInstaller so
end users do not need Python at all (external tools — vina/PyMOL — are still
called by path, so install those via Path B first).

```bash
pip install pyinstaller

pyinstaller --noconfirm --clean --windowed --name DockFlow \
    --collect-all dockflow_core --collect-all dockflow_gui \
    dockflow_gui/app.py
```

Result:

- **Windows**: `dist/DockFlow/DockFlow.exe` (zip the folder to distribute)
- **Linux**: `dist/DockFlow/DockFlow`
- **macOS**: `dist/DockFlow.app` (sign with `codesign` to avoid Gatekeeper
  prompts when distributing)

Notes:

- Data-driven imports (RDKit/Meeko) are picked up by `--collect-all
  dockflow_core`; if RDKit submodules are reported missing add
  `--collect-all rdkit --collect-all meeko`.
- The CLI can be frozen similarly with `dockflow_core/cli.py` (drop
  `--windowed`).
- Keep `vina.exe`/`pymol` available on `PATH`, or set `DOCKFLOW_VINA` /
  `DOCKFLOW_PYMOL` system-wide — frozen apps read the same variables.

## 9. Verification checklist

After any install path, confirm in a terminal:

```bash
dockflow info                       # environment report (see below)
pytest -m "not network and not gui" # offline test suite (133 tests)
python -c "import dockflow_bindings as b; assert callable(b.kabsch_rmsd)"  # accelerator
dockflow-gui                        # GUI launches
```

A healthy full-stack `dockflow info` shows:

```text
rdkit              : 2026.03.x
meeko              : 0.8.x
openbabel          : 3.1.x
vina               : 1.2.x          (or "not installed" - fine if CLI backend present)
backend:python     : 1.2.x
backend:cli        : AutoDock Vina v1.2.7
prep engine        : openbabel
```

Red flags and fixes: see the next section.

## 10. Troubleshooting

### All platforms

| symptom | fix |
|---|---|
| `no Vina backend available` | `pip install vina` (Linux) or install the CLI (§4) — check `dockflow info` |
| `meeko is required` / `ModuleNotFoundError: gemmi` | `pip install ".[prep]"` — gemmi is included since 0.1.1 |
| bindings import fails | rebuild: `pip install ./bindings`; ensure a C++17 compiler is installed; check §7.2 |
| slow analysis on huge pose files | install the bindings (§7.2) — look for `dockflow_bindings` in `dockflow info` |
| renders look schematic | matplotlib fallback is active — install `pymol-open-source` (conda-forge) |

### Linux

| symptom | fix |
|---|---|
| GUI: `libEGL.so.1: cannot open shared object file` | `sudo apt install libegl1 libxkbcommon0 libxkbcommon-x11-0 libgl1` |
| GUI on headless server | `QT_QPA_PLATFORM=offscreen dockflow-gui` (or use the CLI) |
| vina python wheel install fails on non-x86_64 | use the CLI backend: `conda install -c bioconda autodock-vina` |
| `run_dockflow.sh` → "Python found, but DockFlow is not installed in it" | expected on very first run — choose `[1]` (full conda setup) or `[2]` (quick pip setup + Vina download), or run `bash scripts/install_tools.sh` yourself |
| quick setup pip fails with "externally-managed-environment" | the launcher retries automatically with `--break-system-packages`; if that also fails use the full conda setup (`[1]`) |

### Windows

| symptom | fix |
|---|---|
| double-click `run_dockflow.bat` → `ModuleNotFoundError: No module named 'yaml'` | normal on very first run — nothing is installed yet. The launcher offers a guided setup: press `[1]` for the full conda setup or `[2]` for quick pip setup. Or run one of these yourself in the repo folder: `powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1` or `python -m pip install ".[prep,gui,viz]"` |
| `pip install vina` tries to compile for minutes then fails | expected — no Windows wheels; use `vina.exe` (§4). The CLI backend is selected automatically |
| PowerShell blocks `install_tools.ps1` | run `powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1` |
| `cl.exe` not found when building bindings | install VS Build Tools with the C++ workload; run from "x64 Native Tools Command Prompt" |
| path too long errors during pip build | `git config --system core.longpaths true` + enable Windows long paths, or install to `C:\df` |
| SmartScreen warns on the frozen exe | normal for unsigned binaries — "More info → Run anyway", or sign with `signtool` |

### macOS

| symptom | fix |
|---|---|
| `xcode-select: note: command line tools required` | `xcode-select --install` |
| conda env created for wrong architecture | install the Miniforge matching your chip (arm64 vs x86_64); `file $(which python)` shows the arch |
| vina python bindings fail to build | use the bioconda CLI: `conda install -c bioconda autodock-vina`, or let the launcher's quick setup download the official `vina_1.2.7_mac_*` binary |
| "cannot be opened because the developer cannot be verified" (frozen app) | System Settings → Privacy & Security → Open Anyway, or `xattr -dr com.apple.quarantine DockFlow.app` |
| `run_dockflow.command` does nothing on double-click | `chmod +x run_dockflow.command` once (zips strip the executable bit) |
| `run_dockflow.command` → "Python found, but DockFlow is not installed in it" | expected on very first run — choose `[1]` (full conda setup) or `[2]` (quick pip setup + Vina download) |

## 11. What CI does

`.github/workflows/build.yml` runs on every push/PR and on `v*` tags:

| job | what | platforms |
|---|---|---|
| `lint` | `ruff check .` | Linux |
| `test` | offline pytest matrix | Ubuntu, macOS, Windows × Python 3.10 & 3.12 |
| `gui` | PyQt6 offscreen GUI smoke tests | all three OS |
| `bindings` | cibuildwheel builds of the C++ accelerator (cp310–cp312, macOS x86_64+arm64) | all three OS |
| `docker` | builds + smoke-tests the pinned image | Linux |
| `release` | on `v*` tags: sdist + wheel + accelerator wheels attached to the GitHub release | Linux |

So every commit is proven to install, import, test and compile on all three
operating systems.

---

*Related documents: [README.md](README.md) · [USER_GUIDE.md](USER_GUIDE.md) ·
[CHANGELOG.md](CHANGELOG.md)*
