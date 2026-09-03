#!/usr/bin/env bash
# =============================================================================
# run_dockflow.sh - Linux/macOS launcher for DockFlow-Automator.
#
#   ./run_dockflow.sh                GUI (default)
#   ./run_dockflow.sh --cli info     CLI: environment report
#   ./run_dockflow.sh --cli run --config examples/configs/hiv1_protease_example.yaml
#
# Resolution order:
#   1. an activated conda/mamba env that already provides dockflow
#   2. a conda env named "dockflow" in the usual locations
#   3. the system python - and if DockFlow is not installed there yet,
#      this launcher offers a GUIDED SETUP:
#        [1] full setup  (conda env + Vina + PyMOL, ~10 min)
#        [2] quick setup (pip into this Python + vina engine download)
#
# A "vina" executable dropped in this folder is detected and used
# automatically (DOCKFLOW_VINA).
#
# Full setup (can also be run directly, once):
#   bash scripts/install_tools.sh
# =============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV_NAME="${DOCKFLOW_ENV:-dockflow}"

blue()  { printf '\033[1;34m[dockflow]\033[0m %s\n' "$*"; }
err()   { printf '\033[1;31m[dockflow] error:\033[0m %s\n' "$*" >&2; }

# ---------------------------------------------------------------- environment
find_env_python() {
    # already inside an activated environment?
    if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/dockflow-gui" ]]; then
        echo "${CONDA_PREFIX}"
        return 0
    fi
    for base in \
        "${HOME}/miniforge3" "${HOME}/mambaforge" "${HOME}/miniconda3" \
        "${HOME}/anaconda3" "/opt/miniforge3" "/opt/conda" "/usr/local/miniforge3" \
        "/opt/homebrew/Caskroom/miniforge/base" "/usr/local/Caskroom/miniforge/base"
    do
        if [[ -x "$base/envs/$ENV_NAME/bin/dockflow-gui" ]]; then
            echo "$base/envs/$ENV_NAME"
            return 0
        fi
    done
    return 1
}

if ENVDIR="$(find_env_python)"; then
    blue "using environment: $ENVDIR"
    if [[ "${1:-}" == "--cli" ]]; then
        shift
        exec "$ENVDIR/bin/dockflow" "$@"
    fi
    exec "$ENVDIR/bin/dockflow-gui" "$@"
fi

# ---------------------------------------------------------- system python path
blue "no conda environment \"$ENV_NAME\" found - trying system python"
PY="$(command -v python3 || true)"
if [[ -z "$PY" ]]; then
    err "no python3 on PATH."
    echo "  Run the full setup once - it installs its own Python:" >&2
    echo "      bash scripts/install_tools.sh" >&2
    echo "  (needs Miniforge: https://github.com/conda-forge/miniforge/releases)" >&2
    exit 1
fi

# a "vina" executable dropped next to this launcher is picked up automatically
if [[ -x "$ROOT/vina" ]]; then
    export DOCKFLOW_VINA="$ROOT/vina"
fi

# ------------------------------------------------------------ dependency check
if [[ "${1:-}" == "--cli" ]]; then
    "$PY" -c "import yaml, numpy, requests" >/dev/null 2>&1 && DEPS_OK=1 || DEPS_OK=0
else
    "$PY" -c "import yaml, numpy, requests, PyQt6" >/dev/null 2>&1 && DEPS_OK=1 || DEPS_OK=0
fi
if [[ "$DEPS_OK" == "1" ]]; then
    if [[ "${1:-}" == "--cli" ]]; then
        shift
        exec "$PY" -m dockflow_core.cli "$@"
    fi
    if ! command -v vina >/dev/null 2>&1 && [[ -z "${DOCKFLOW_VINA:-}" ]]; then
        blue "note: no Vina engine detected yet - the GUI will start, but to"
        blue "      dock, run the setup: bash scripts/install_tools.sh"
    fi
    exec "$PY" -m dockflow_gui "$@"
fi

# ---------------------------------------------------------------- guided setup
blue "Python found ($PY), but DockFlow is not installed in it."
echo
echo "  [1] Full setup  - conda environment + Vina + PyMOL - recommended, ~10 min"
echo "  [2] Quick setup - pip install into this Python + vina download, ~2 min"
echo "  [3] Exit"
echo
choice=3
read -r -p "Choose an option [1-3]: " choice || choice=3
case "$choice" in
    1)
        blue "running the full setup: scripts/install_tools.sh"
        if ! bash scripts/install_tools.sh; then
            err "the full setup failed - see the messages above."
            echo "  Trying the quick pip setup instead..." >&2
            choice=2
        else
            blue "full setup finished - restarting the app..."
            exec "$0" "$@"
        fi
        ;;
esac
if [[ "$choice" != "2" ]]; then
    exit 0
fi

# ----------------------------------------------------------------- quick setup
blue "quick setup: pip install \".[prep,gui,viz]\""
if ! "$PY" -m pip install ".[prep,gui,viz]"; then
    blue "plain pip install failed - retrying with --break-system-packages"
    blue "(PEP 668 'externally managed environment', e.g. Homebrew python)..."
    if ! "$PY" -m pip install --break-system-packages ".[prep,gui,viz]"; then
        err "pip install failed - check the messages above."
        echo "  Recommended instead: bash scripts/install_tools.sh (full setup)" >&2
        exit 1
    fi
fi

# fetch the official Vina engine so docking works out of the box
if [[ -z "${DOCKFLOW_VINA:-}" && ! -x "$ROOT/vina" ]]; then
    os="$(uname -s)"
    arch="$(uname -m)"
    case "$os/$arch" in
        Darwin/arm64)   asset="vina_1.2.7_mac_aarch64" ;;
        Darwin/x86_64)  asset="vina_1.2.7_mac_x86_64" ;;
        Darwin/aarch64) asset="vina_1.2.7_mac_aarch64" ;;
        Linux/aarch64)  asset="vina_1.2.7_linux_aarch64" ;;
        *)              asset="vina_1.2.7_linux_x86_64" ;;
    esac
    blue "downloading the AutoDock Vina 1.2.7 engine - $asset..."
    url="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/$asset"
    if (curl -fsSL --max-time 120 "$url" -o "$ROOT/vina" 2>/dev/null \
            || wget -q --timeout=120 "$url" -O "$ROOT/vina" 2>/dev/null); then
        chmod +x "$ROOT/vina"
        blue "vina saved in this folder - docking is ready."
    else
        rm -f "$ROOT/vina"
        blue "note: vina could not be downloaded. The app will still start,"
        blue "      but to dock, run the full setup (option [1]) or download"
        blue "      $asset from https://github.com/ccsb-scripps/AutoDock-Vina/releases"
        blue "      rename it to 'vina', put it in this folder, restart the app."
    fi
fi
if [[ -x "$ROOT/vina" ]]; then
    export DOCKFLOW_VINA="$ROOT/vina"
fi

blue "quick setup finished - starting the app..."
if [[ "${1:-}" == "--cli" ]]; then
    shift
    exec "$PY" -m dockflow_core.cli "$@"
fi
exec "$PY" -m dockflow_gui "$@"
