#!/usr/bin/env bash
# DockFlow-Automator container entrypoint.
#
# All arguments are forwarded to the `dockflow` CLI.  The special first
# argument "gui" launches the PyQt6 desktop application (requires X11
# forwarding: -e DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix), and "shell"
# starts an interactive bash session inside the environment.
#
# NOTE: micromamba environments live in ${MAMBA_ROOT_PREFIX}/envs/<name>
# (default /opt/conda/envs/dockflow) - "env", not "envs", is the classic
# path typo that breaks every CLI call in the container.
set -euo pipefail

ENV_BIN="/opt/conda/envs/dockflow/bin"
export PATH="${ENV_BIN}:${PATH}"
export DOCKFLOW_HOME="${DOCKFLOW_HOME:-/data}"

case "${1:-}" in
    gui)
        shift
        exec "${ENV_BIN}/dockflow-gui" "$@"
        ;;
    shell)
        # pass through any args, e.g.: shell -c "command -v vina"
        exec /bin/bash "$@"
        ;;
    dockflow)
        # allow the explicit "dockflow ..." form too
        shift
        exec "${ENV_BIN}/dockflow" "$@"
        ;;
    *)
        exec "${ENV_BIN}/dockflow" "$@"
        ;;
esac
