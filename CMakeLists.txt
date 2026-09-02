# =============================================================================
# DockFlow-Automator - top-level CMake build.
#
# This "superbuild" compiles the optional C++ accelerator bindings
# (bindings/dockflow_bindings) that provide fast PDBQT parsing, grid-box
# computation, Kabsch RMSD and contact evaluation to Python.
#
# The Python application itself installs through pip:
#     pip install .
#     pip install ./bindings     # scikit-build-core wheel for the bindings
#
# Standalone CMake usage (drops the module next to the sources):
#     cmake -B build -DPython3_EXECUTABLE=$(which python)
#     cmake --build build
#     cmake --install build --prefix <site-packages dir>
#
# OpenBabel can also be built from source here for full reproducibility
# (OFF by default; the Python wheels or conda-forge are usually easier):
#     cmake -B build -DDOCKFLOW_BUILD_OPENBABEL=ON
# =============================================================================
cmake_minimum_required(VERSION 3.18)

project(DockFlowAutomator
    VERSION 0.1.0
    DESCRIPTION "Unified, automated molecular docking platform"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(DOCKFLOW_BUILD_BINDINGS "Build the pybind11 accelerator module" ON)
option(DOCKFLOW_BUILD_OPENBABEL "Also build OpenBabel from source (heavy)" OFF)

# ---------------------------------------------------------------------------
# Accelerator bindings (pybind11)
# ---------------------------------------------------------------------------
if(DOCKFLOW_BUILD_BINDINGS)
    add_subdirectory(bindings)
endif()

# ---------------------------------------------------------------------------
# Optional OpenBabel superbuild (ExternalProject)
# ---------------------------------------------------------------------------
if(DOCKFLOW_BUILD_OPENBABEL)
    include(ExternalProject)
    ExternalProject_Add(openbabel
        GIT_REPOSITORY https://github.com/openbabel/openbabel.git
        GIT_TAG openbabel-3-1-1
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/openbabel-install
            -DBUILD_GUI=OFF
            -DPython_BINDINGS=ON
            -DENABLE_TESTS=OFF
        INSTALL_DIR ${CMAKE_BINARY_DIR}/openbabel-install)
    message(STATUS "OpenBabel will be built into ${CMAKE_BINARY_DIR}/openbabel-install")
    message(STATUS "Point BABEL_LIBDIR / LD_LIBRARY_PATH there at runtime.")
endif()

message(STATUS "")
message(STATUS "DockFlow-Automator ${PROJECT_VERSION} configuration:")
message(STATUS "  C++ compiler   : ${CMAKE_CXX_COMPILER}")
message(STATUS "  Bindings       : ${DOCKFLOW_BUILD_BINDINGS}")
message(STATUS "  OpenBabel superbuild : ${DOCKFLOW_BUILD_OPENBABEL}")
message(STATUS "")
