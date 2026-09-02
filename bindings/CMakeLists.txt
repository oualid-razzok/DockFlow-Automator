# =============================================================================
# dockflow_bindings - pybind11 accelerator module for DockFlow-Automator.
#
# Built either through scikit-build-core (pip install ./bindings) or through
# a plain CMake invocation:
#     cmake -B build -DPython3_EXECUTABLE=$(which python)
#     cmake --build build
# =============================================================================
cmake_minimum_required(VERSION 3.18)

project(dockflow_bindings
    VERSION 0.1.0
    DESCRIPTION "C++ accelerator kernels for DockFlow-Automator"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(DOCKFLOW_BINDINGS_WERROR "Treat warnings as errors" OFF)
if(DOCKFLOW_BINDINGS_WERROR)
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# ---------------------------------------------------------------------------
# pybind11: use the installed package when available, otherwise fetch it.
# ---------------------------------------------------------------------------
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.12.0
    )
    FetchContent_MakeAvailable(pybind11)
    message(STATUS "pybind11 fetched from upstream (v2.12.0)")
else()
    message(STATUS "pybind11 found via CONFIG: ${pybind11_DIR}")
endif()

pybind11_add_module(dockflow_bindings MODULE
    src/dockflow_bindings.cpp
)

# scikit-build-core installs into the wheel root; a plain CMake install goes
# to the site-packages of the interpreter used for configuration.
if(SKBUILD)
    install(TARGETS dockflow_bindings LIBRARY DESTINATION .)
else()
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c
        "import sysconfig, os; print(sysconfig.get_paths()['purelib'])"
        OUTPUT_VARIABLE DOCKFLOW_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE DOCKFLOW_SITE_RESULT
    )
    if(DOCKFLOW_SITE_RESULT EQUAL 0 AND DOCKFLOW_SITE_PACKAGES)
        install(TARGETS dockflow_bindings LIBRARY DESTINATION ${DOCKFLOW_SITE_PACKAGES})
    else()
        install(TARGETS dockflow_bindings LIBRARY DESTINATION lib)
    endif()
endif()

message(STATUS "dockflow_bindings ${PROJECT_VERSION} configured for "
        "Python ${Python3_VERSION}")
