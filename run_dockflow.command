#!/usr/bin/env bash
# =============================================================================
# run_dockflow.command - macOS launcher for DockFlow-Automator.
#
# Double-click this file in Finder to start the desktop GUI (it opens a
# Terminal window), or invoke it from a shell like a normal script:
#
#   ./run_dockflow.command                 GUI (default)
#   ./run_dockflow.command --cli info      CLI: environment report
#   ./run_dockflow.command --cli run --config examples/configs/hiv1_protease_example.yaml
#
# Resolution order:
#   1. an activated conda/mamba env that already provides dockflow
#   2. a conda env named "dockflow" in the usual macOS locations
#   3. the Homebrew / system python (dockflow-automator pip-installed)
#
# NOTE: files downloaded as part of a .zip may lose the executable bit.
# If double-clicking fails, run once in Terminal:
#   chmod +x run_dockflow.command
# =============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV_NAME="${DOCKFLOW_ENV:-dockflow}"

find_env_python() {
    # already inside an activated environment?
    if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/dockflow-gui" ]]; then
        echo "${CONDA_PREFIX}"
        return 0
    fi
    for base in \
        "${HOME}/miniforge3" "${HOME}/mambaforge" "${HOME}/miniconda3" \
        "${HOME}/anaconda3" "/opt/homebrew/Caskroom/miniforge/base" \
        "/usr/local/Caskroom/miniforge/base" "/opt/miniforge3"
    do
        if [[ -x "$base/envs/$ENV_NAME/bin/dockflow-gui" ]]; then
            echo "$base/envs/$ENV_NAME"
            return 0
        fi
    done
    return 1
}

if ENVDIR="$(find_env_python)"; then
    printf '\033[1;34m[dockflow]\033[0m using environment: %s\n' "$ENVDIR"
    if [[ "${1:-}" == "--cli" ]]; then
        shift
        exec "$ENVDIR/bin/dockflow" "$@"
    else
        exec "$ENVDIR/bin/dockflow-gui" "$@"
    fi
fi

printf '\033[1;34m[dockflow]\033[0m no conda environment "%s" found - trying system python\n' "$ENV_NAME"
PY="$(command -v python3 || true)"
if [[ -z "$PY" ]]; then
    printf '\033[1;31m[dockflow] error:\033[0m no python3 on PATH.\n' >&2
    echo "  Install Miniforge and run: bash scripts/install_tools.sh" >&2
    read -r -p "Press return to close..." _
    exit 1
fi
if [[ "${1:-}" == "--cli" ]]; then
    shift
    exec "$PY" -m dockflow_core.cli "$@"
fi
exec "$PY" -m dockflow_gui "$@"
