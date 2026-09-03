#!/usr/bin/env bash
# =============================================================================
# run_dockflow.sh - Linux launcher for DockFlow-Automator.
#
#   ./run_dockflow.sh                GUI (default)
#   ./run_dockflow.sh --cli info     CLI: environment report
#   ./run_dockflow.sh --cli run --config examples/configs/hiv1_protease_example.yaml
#
# Resolution order:
#   1. an activated conda/mamba env that already provides dockflow
#   2. a conda env named "dockflow" in the usual locations
#   3. the system python (dockflow-automator must be pip-installed)
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
        "${HOME}/anaconda3" "/opt/miniforge3" "/opt/conda" "/usr/local/miniforge3"
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
if ! command -v python3 >/dev/null 2>&1; then
    printf '\033[1;31m[dockflow] error:\033[0m no python3 on PATH.\n' >&2
    echo "  Run: bash scripts/install_tools.sh   (creates the environment)" >&2
    exit 1
fi
if [[ "${1:-}" == "--cli" ]]; then
    shift
    exec python3 -m dockflow_core.cli "$@"
fi
exec python3 -m dockflow_gui "$@"
