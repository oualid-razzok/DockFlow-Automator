# =============================================================================
# DockFlow-Automator reproducible environment.
#
# Guarantees the heavy scientific stack that pip cannot always provide:
# OpenBabel (python bindings), open-source PyMOL (headless rendering),
# the AutoDock Vina executable and swig (to build the vina python bindings).
#
# Build:  docker build -f docker/Dockerfile -t dockflow-automator .
# Run:    docker run --rm -v "$PWD:/data" dockflow-automator info
#         docker run --rm -v "$PWD:/data" dockflow-automator run --config /data/run.yaml
# GUI:    docker run --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
#             -v "$PWD:/data" dockflow-automator gui   (X11 forwarding)
# =============================================================================
FROM mambaorg/micromamba:jammy-1.5.10 AS base

# System libraries: OpenGL/rendering for headless PyMOL ray tracing,
# X11 client libs for optional GUI forwarding, compilers for the bindings.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libgl1 \
        libglu1-mesa \
        libxrender1 \
        libxext6 \
        libxi6 \
        xvfb \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Conda environment with the scientific stack
# ---------------------------------------------------------------------------
COPY --chown=${MAMBA_USER}:${MAMBA_USER} docker/env.yaml /tmp/env.yaml
RUN micromamba create -y -n dockflow -c conda-forge -c bioconda -f /tmp/env.yaml \
    && micromamba clean --all -y \
    && rm /tmp/env.yaml

# ---------------------------------------------------------------------------
# Application + accelerator bindings
# ---------------------------------------------------------------------------
COPY --chown=${MAMBA_USER}:${MAMBA_USER} . /opt/dockflow
RUN micromamba run -n dockflow pip install --no-cache-dir \
        "/opt/dockflow[prep,viz]" \
    && micromamba run -n dockflow pip install --no-cache-dir \
        "/opt/dockflow/bindings" \
    && micromamba run -n dockflow python -c "import dockflow_core, dockflow_bindings" \
    || echo "bindings optional - core still usable"

# ---------------------------------------------------------------------------
# Runtime image
# ---------------------------------------------------------------------------
FROM base
ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    DOCKFLOW_HOME=/data \
    MPLBACKEND=Agg \
    QT_QPA_PLATFORM=offscreen
WORKDIR /data
VOLUME ["/data"]

COPY --from=base --chown=${MAMBA_USER}:${MAMBA_USER} /opt/dockflow/docker/entrypoint.sh \
     /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh \
    && mkdir -p /data \
    && chown -R ${MAMBA_USER} /data /opt/dockflow

USER $MAMBA_USER
LABEL org.opencontainers.image.title="DockFlow-Automator" \
      org.opencontainers.image.description="Unified, automated molecular docking platform" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/oualid-razzok/DockFlow-Automator"

HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD micromamba run -n dockflow dockflow --version || exit 1

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["info"]
