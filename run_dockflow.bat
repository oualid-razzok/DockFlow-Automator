@echo off
REM =============================================================================
REM run_dockflow.bat - Windows launcher for DockFlow-Automator.
REM
REM Double-click to start the desktop GUI, or run from a terminal:
REM    run_dockflow.bat                     GUI (default)
REM    run_dockflow.bat --cli info          CLI: environment report
REM    run_dockflow.bat --cli run --config examples\configs\hiv1_protease_example.yaml
REM
REM Resolution order:
REM    1. an activated conda environment (CONDA_PREFIX with dockflow installed)
REM    2. a conda environment named "dockflow" in the usual install locations
REM    3. the system Python (dockflow-automator must be pip-installed)
REM
REM Create the environment with:  powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
REM =============================================================================
setlocal
cd /d "%~dp0"

set "DOCKFLOW_ENV_NAME=dockflow"
set "ENVDIR="

REM -- 1. already inside an activated environment? -----------------------------
if defined CONDA_PREFIX (
    if exist "%CONDA_PREFIX%\Scripts\dockflow-gui.exe" set "ENVDIR=%CONDA_PREFIX%"
)

REM -- 2. look for a conda env named "dockflow" --------------------------------
if not defined ENVDIR (
    for %%B in (
        "%USERPROFILE%\miniforge3"
        "%USERPROFILE%\mambaforge"
        "%USERPROFILE%\miniconda3"
        "%USERPROFILE%\anaconda3"
        "C:\ProgramData\miniforge3"
        "C:\ProgramData\miniconda3"
        "C:\ProgramData\anaconda3"
    ) do (
        if not defined ENVDIR (
            if exist "%%~fB\envs\%DOCKFLOW_ENV_NAME%\Scripts\dockflow-gui.exe" (
                set "ENVDIR=%%~fB\envs\%DOCKFLOW_ENV_NAME%"
            )
        )
    )
)

REM -- dispatch ------------------------------------------------------------------
if "%~1"=="--cli" goto :cli
goto :gui

REM -- CLI mode -------------------------------------------------------------------
:cli
if defined ENVDIR (
    echo [dockflow] using environment: %ENVDIR%
    "%ENVDIR%\Scripts\dockflow.exe" %2 %3 %4 %5 %6 %7 %8 %9
    goto :end
)
echo [dockflow] no conda environment "%DOCKFLOW_ENV_NAME%" found - trying system python
where python >nul 2>nul
if errorlevel 1 (
    echo [dockflow] error: no python found on PATH.
    echo    Run scripts\install_tools.ps1 first to create the environment.
    pause
    exit /b 1
)
python -m dockflow_core.cli %2 %3 %4 %5 %6 %7 %8 %9
goto :end

REM -- GUI mode (default) ----------------------------------------------------------
:gui
if defined ENVDIR (
    echo [dockflow] using environment: %ENVDIR%
    "%ENVDIR%\Scripts\dockflow-gui.exe" %*
    goto :end
)
echo [dockflow] no conda environment "%DOCKFLOW_ENV_NAME%" found - trying system python
where python >nul 2>nul
if errorlevel 1 (
    echo [dockflow] error: no python found on PATH.
    echo    Run scripts\install_tools.ps1 first to create the environment.
    pause
    exit /b 1
)
python -m dockflow_gui %*

:end
if errorlevel 1 (
    echo [dockflow] the command exited with an error.
    pause
)
endlocal
