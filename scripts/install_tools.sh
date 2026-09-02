#!/usr/bin/env bash
# =============================================================================
# install_tools.sh - install the heavy external tools DockFlow orchestrates.
#
# Creates a conda/mamba environment "dockflow" with:
#   - openbabel            (format conversion, receptor preparation)
#   - pymol-open-source    (3D visualization / rendering)
#   - autodock-vina        (Vina CLI executable, bioconda)
# and then pip-installs this repository + the vina python bindings (needs
# swig, provided by conda) + the C++ accelerator bindings.
#
# Usage:
#   bash scripts/install_tools.sh            # use micromamba if found
#   bash scripts/install_tools.sh --conda    # force conda
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${DOCKFLOW_ENV:-dockflow}"
PYTHON_VERSION="${DOCKFLOW_PYTHON:-3.10}"

say() { printf '\033[1;34m[dockflow]\033[0m %s\n' "$*"; }
die() { printf '\033[1;31m[dockflow] error:\033[0m %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# pick a package manager
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "--conda" ]]; then
    CONDA_BIN="$(command -v conda || true)"
else
    CONDA_BIN="$(command -v micromamba || command -v mamba || command -v conda || true)"
fi
[[ -n "$CONDA_BIN" ]] || die "no conda/mamba/micromamba found; install miniforge first:
  https://github.com/conda-forge/miniforge"

say "using package manager: $CONDA_BIN"
case "$CONDA_BIN" in
    *micromamba*) CREATE=(micromamba create -y -n "$ENV_NAME" -c conda-forge -c bioconda) ; RUN=(micromamba run -n "$ENV_NAME") ;;
    *mamba*)      CREATE=(mamba create -y -n "$ENV_NAME" -c conda-forge -c bioconda) ;      RUN=(conda run -n "$ENV_NAME") ;;
    *)            CREATE=(conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda) ;      RUN=(conda run -n "$ENV_NAME") ;;
esac

# ---------------------------------------------------------------------------
# environment
# ---------------------------------------------------------------------------
say "creating environment '$ENV_NAME' (python $PYTHON_VERSION, openbabel, pymol, vina)"
"${CREATE[@]}" \
    "python=$PYTHON_VERSION" \
    openbabel \
    pymol-open-source \
    autodock-vina \
    swig \
    pip

say "environment ready; installing DockFlow packages"
"${RUN[@]}" pip install --quiet "${ROOT}[prep,engine,gui,viz]"
if [[ -d "$ROOT/bindings" ]]; then
    say "building the C++ accelerator bindings (needs a C++17 compiler)"
    "${RUN[@]}" pip install --quiet "$ROOT/bindings" || \
        say "bindings build failed (optional) - continuing without acceleration"
fi

say "verifying installation"
"${RUN[@]}" python -c "import openbabel; print('openbabel OK')" || true
"${RUN[@]}" python -c "import pymol; print('pymol OK')" || true
"${RUN[@]}" vina --version || true

cat <<'EOF'

Done. Activate and start:

  conda activate dockflow          # (or: micromamba activate dockflow)
  dockflow info                    # environment report
  dockflow gui                     # desktop application
  dockflow run --config my_run.yaml

EOF
