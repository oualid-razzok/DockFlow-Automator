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
REM    3. the system Python - and if DockFlow is not installed there yet,
REM       this launcher offers a GUIDED SETUP:
REM         [1] full setup  (conda env + Vina + PyMOL, ~10 min)
REM         [2] quick setup (pip into this Python + auto-download of vina.exe)
REM
REM A vina.exe dropped in this folder is detected and used automatically.
REM
REM Full setup (can also be run directly, once):
REM    powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
REM =============================================================================
setlocal
cd /d "%~dp0"

set "DOCKFLOW_ENV_NAME=dockflow"
set "ENVDIR="
set "MODE=gui"
if /i "%~1"=="--cli" set "MODE=cli"

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

REM -- launch from a conda environment ------------------------------------------
if defined ENVDIR (
    echo [dockflow] using environment: %ENVDIR%
    if "%MODE%"=="cli" (
        "%ENVDIR%\Scripts\dockflow.exe" %2 %3 %4 %5 %6 %7 %8 %9
    ) else (
        "%ENVDIR%\Scripts\dockflow-gui.exe" %*
    )
    goto :end
)

REM -- 3. system python fallback -------------------------------------------------
set "PY="
where python >nul 2>nul
if not errorlevel 1 set "PY=python"
if not defined PY (
    where py >nul 2>nul
    if not errorlevel 1 set "PY=py -3"
)
if not defined PY (
    echo [dockflow] no usable Python found on this computer.
    echo.
    echo    Either run the full setup once - it installs its own Python:
    echo        powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
    echo    needs Miniforge: https://github.com/conda-forge/miniforge/releases
    echo.
    echo    or install Python from https://www.python.org/downloads/
    echo    then double-click this file again.
    pause
    exit /b 1
)

REM guard against the Windows Store python stub - exists but cannot run
%PY% -c "import sys" >nul 2>nul
if errorlevel 1 (
    echo [dockflow] "%PY%" exists but cannot run scripts - Windows Store alias?
    echo    Run the full setup instead - it installs its own Python:
    echo        powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
    echo    or install real Python from https://www.python.org/downloads/
    pause
    exit /b 1
)

REM a vina.exe dropped next to this launcher is picked up automatically
if exist "%~dp0vina.exe" set "DOCKFLOW_VINA=%~dp0vina.exe"

REM -- dependency check ----------------------------------------------------------
if "%MODE%"=="cli" (
    %PY% -c "import yaml, numpy, requests" >nul 2>nul
) else (
    %PY% -c "import yaml, numpy, requests, PyQt6" >nul 2>nul
)
if not errorlevel 1 goto %MODE%

REM -- guided setup ---------------------------------------------------------------
echo [dockflow] Python found, but DockFlow is not installed in it.
echo.
echo   [1] Full setup  - conda environment + Vina + PyMOL - recommended, ~10 min
echo   [2] Quick setup - pip install into this Python + vina.exe download, ~2 min
echo   [3] Exit
choice /c 123 /n /m "Choose an option [1-3]: "
if errorlevel 3 goto :end
if errorlevel 2 goto :quicksetup

echo [dockflow] running the full setup: scripts\install_tools.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File "scripts\install_tools.ps1"
if errorlevel 1 (
    echo [dockflow] the full setup failed - see the messages above.
    echo [dockflow] trying the quick pip setup instead...
    goto :quicksetup
)
echo [dockflow] full setup finished - restarting the app...
endlocal
"%~f0" %*
exit /b 0

:quicksetup
echo [dockflow] quick setup: pip install ".[prep,gui,viz]"
%PY% -m pip install ".[prep,gui,viz]"
if errorlevel 1 (
    echo [dockflow] pip install failed - check the messages above.
    pause
    exit /b 1
)

REM fetch the official Vina engine so docking works out of the box
if not exist "%~dp0vina.exe" (
    echo [dockflow] downloading the AutoDock Vina 1.2.7 engine - vina.exe...
    powershell -NoProfile -Command "try { [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri 'https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe' -OutFile 'vina.exe' -UseBasicParsing } catch { exit 1 }" >nul 2>nul
    if exist "%~dp0vina.exe" (
        echo [dockflow] vina.exe saved in this folder - docking is ready.
    ) else (
        echo [dockflow] note: vina.exe could not be downloaded. The app will
        echo    still start, but to dock you need it - download manually from
        echo    https://github.com/ccsb-scripps/AutoDock-Vina/releases
        echo    rename to vina.exe, put it in this folder, restart the app.
    )
)
if exist "%~dp0vina.exe" set "DOCKFLOW_VINA=%~dp0vina.exe"
goto %MODE%

:cli
%PY% -m dockflow_core.cli %2 %3 %4 %5 %6 %7 %8 %9
goto :end

:gui
where vina.exe >nul 2>nul
if errorlevel 1 (
    if not exist "%~dp0vina.exe" (
        echo [dockflow] note: no Vina engine detected yet - the GUI will start,
        echo    but to dock, drop vina.exe into this folder or run the full
        echo    setup: powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
    )
)
%PY% -m dockflow_gui %*
goto :end

:end
if errorlevel 1 (
    echo [dockflow] the command exited with an error.
    pause
)
endlocal
