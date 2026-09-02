#!/usr/bin/env bash
# DockFlow-Automator container entrypoint.
#
# All arguments are forwarded to the `dockflow` CLI.  The special first
# argument "gui" launches the PyQt6 desktop application (requires X11
# forwarding: -e DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix), and "shell"
# starts an interactive bash session inside the environment.
set -euo pipefail

export PATH="/opt/conda/env/dockflow/bin:${PATH}"
export DOCKFLOW_HOME="${DOCKFLOW_HOME:-/data}"

case "${1:-}" in
    gui)
        shift
        exec /opt/conda/env/dockflow/bin/dockflow-gui "$@"
        ;;
    shell)
        exec /bin/bash
        ;;
    *)
        exec /opt/conda/env/dockflow/bin/dockflow "$@"
        ;;
esac
