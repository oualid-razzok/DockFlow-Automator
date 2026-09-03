# =============================================================================
# install_tools.ps1 - Windows bootstrap for DockFlow-Automator.
#
# Creates a conda environment "dockflow" with:
#   - openbabel          (format conversion, receptor prep; conda-forge)
#   - pymol-open-source  (3D visualization; conda-forge)
#   - the official AutoDock Vina 1.2.7 Windows executable (renamed vina.exe,
#     placed in the environment so `dockflow info` finds it on PATH)
# and then pip-installs the repository with the [prep,gui,viz] extras plus the
# optional C++ accelerator bindings (built with MSVC if available).
#
# NOTE: the vina *python bindings* (pip package "vina") have no Windows
# wheels; docking on Windows uses the equally capable CLI backend, which the
# engine auto-detects. No functionality is lost: backends produce identical
# PDBQT/CSV outputs.
#
# Usage (PowerShell):
#   powershell -ExecutionPolicy Bypass -File scripts\install_tools.ps1
#
# Optional flags:
#   -EnvName dockflow        environment name
#   -PythonVersion 3.10      python in the env
#   -SkipVinaDownload        skip fetching vina.exe (you already have it)
# =============================================================================
param(
    [string]$EnvName = "dockflow",
    [string]$PythonVersion = "3.10",
    [switch]$SkipVinaDownload
)

# NOTE: 'Continue' (not 'Stop') because conda/pip/vina write progress to
# stderr, and under $ErrorActionPreference='Stop' PowerShell 5.1 turns any
# redirected native stderr line into a terminating error. Failures are
# instead handled through explicit $LASTEXITCODE checks below.
$ErrorActionPreference = "Continue"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$Root = Split-Path -Parent $PSScriptRoot   # repository root (scripts\..)
$VinaVersion = "1.2.7"

function Write-Step { param([string]$Message)
    Write-Host "[dockflow] $Message" -ForegroundColor Cyan }
function Write-Ok    { param([string]$Message)
    Write-Host "[dockflow] $Message" -ForegroundColor Green }
function Write-Warn2 { param([string]$Message)
    Write-Host "[dockflow] $Message" -ForegroundColor Yellow }
function Die { param([string]$Message)
    Write-Host "[dockflow] error: $Message" -ForegroundColor Red
    exit 1 }

# -----------------------------------------------------------------------------
# 1. locate a conda-family package manager (miniforge recommended)
# -----------------------------------------------------------------------------
$CondaCandidates = @()
if ($env:CONDA_EXE)     { $CondaCandidates += $env:CONDA_EXE }
if ($env:MAMBA_EXE)     { $CondaCandidates += $env:MAMBA_EXE }
$CondaCandidates += @(
    "$env:USERPROFILE\miniforge3\Scripts\conda.exe",
    "$env:USERPROFILE\miniforge3\condabin\conda.bat",
    "$env:USERPROFILE\mambaforge\Scripts\conda.exe",
    "$env:USERPROFILE\miniconda3\Scripts\conda.exe",
    "$env:USERPROFILE\anaconda3\Scripts\conda.exe",
    "C:\ProgramData\miniforge3\Scripts\conda.exe",
    "C:\ProgramData\miniconda3\Scripts\conda.exe"
)
# also try whatever is on PATH
try { $onPath = Get-Command conda -ErrorAction Stop } catch { $onPath = $null }
if ($onPath) { $CondaCandidates += $onPath.Source }

$CondaExe = $null
foreach ($candidate in $CondaCandidates) {
    if ($candidate -and (Test-Path $candidate)) { $CondaExe = $candidate; break }
}
if (-not $CondaExe) {
    Die "No conda/mamba found. Install Miniforge first:
  https://github.com/conda-forge/miniforge/releases
  (pick 'Miniforge3-Windows-x86_64.exe', then re-run this script)"
}
Write-Step "using conda: $CondaExe"

# -----------------------------------------------------------------------------
# 2. create the environment (conda-forge; bioconda has no Windows channel)
# -----------------------------------------------------------------------------
# base installation: from `conda info --json` (root_prefix), falling back to
# the executable's grand-parent directory (Scripts\conda.exe -> base)
$CondaBase = $null
try {
    $CondaBase = (& $CondaExe info --json 2>$null | ConvertFrom-Json).root_prefix
} catch { $CondaBase = $null }
if (-not $CondaBase) { $CondaBase = Split-Path -Parent (Split-Path -Parent $CondaExe) }
$EnvPath = Join-Path $CondaBase "envs\$EnvName"

Write-Step "creating environment '$EnvName' (python $PythonVersion, openbabel, pymol)"
& $CondaExe create -y -n $EnvName -c conda-forge `
    "python=$PythonVersion" openbabel pymol-open-source pip
if ($LASTEXITCODE -ne 0) { Die "conda create failed (see log above)" }

$EnvPython = Join-Path $EnvPath "python.exe"
$EnvScripts = Join-Path $EnvPath "Scripts"
if (-not (Test-Path $EnvPython)) { Die "environment python not found at $EnvPython" }

# run python inside the env; prints combined stdout/stderr to the console
# and returns $true on success (output goes through Out-Host so it never
# pollutes the function's return value)
function Invoke-EnvPython { param([string[]]$Arguments)
    & $EnvPython @Arguments 2>&1 | ForEach-Object { Write-Host $_ }
    return ($LASTEXITCODE -eq 0)
}

# -----------------------------------------------------------------------------
# 3. official Vina Windows executable -> env\Scripts\vina.exe
#    (the docking engine auto-selects the CLI backend when it finds vina.exe)
# -----------------------------------------------------------------------------
if (-not $SkipVinaDownload) {
    $VinaUrl  = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v$VinaVersion/vina_$($VinaVersion -replace '\.', '_')_win.exe"
    $VinaDest = Join-Path $EnvScripts "vina.exe"
    if (Test-Path $VinaDest) {
        Write-Step "vina.exe already present ($VinaDest)"
    } else {
        Write-Step "downloading AutoDock Vina $VinaVersion for Windows"
        try {
            Invoke-WebRequest -Uri $VinaUrl -OutFile $VinaDest -UseBasicParsing
            Write-Ok "vina.exe -> $VinaDest"
        } catch {
            Write-Warn2 "could not download vina.exe ($($_.Exception.Message))"
            Write-Warn2 "download it manually from https://github.com/ccsb-scripps/AutoDock-Vina/releases"
            Write-Warn2 "rename to vina.exe, place it in $EnvScripts and re-run with -SkipVinaDownload"
        }
    }
} else {
    Write-Step "skipping vina download (-SkipVinaDownload)"
}

# -----------------------------------------------------------------------------
# 4. install the python package (+ extras)
#    engine extra intentionally omitted: no windows wheels for `vina`; the
#    CLI backend (vina.exe above) provides docking with identical results.
# -----------------------------------------------------------------------------
Write-Step "pip-installing DockFlow-Automator [prep,gui,viz]"
$ok = Invoke-EnvPython @("-m", "pip", "install", "--quiet", "$Root[prep,gui,viz]")
if (-not $ok) { Die "pip install of the main package failed" }

# -----------------------------------------------------------------------------
# 5. optional: C++ accelerator bindings (needs MSVC build tools)
# -----------------------------------------------------------------------------
if (Test-Path (Join-Path $Root "bindings")) {
    Write-Step "building the C++ accelerator bindings (requires Visual Studio C++ tools)"
    $ok = Invoke-EnvPython @("-m", "pip", "install", "--quiet", (Join-Path $Root "bindings"))
    if (-not $ok) {
        Write-Warn2 "bindings build failed (optional) - continuing without acceleration"
        Write-Warn2 "install 'Visual Studio Build Tools' (C++ workload) to enable it"
    }
}

# -----------------------------------------------------------------------------
# 6. verify
# -----------------------------------------------------------------------------
Write-Step "verifying the installation"
[void](Invoke-EnvPython @("-c", "import openbabel, openbabel.pybel; print('openbabel OK')"))
[void](Invoke-EnvPython @("-c", "import pymol2; print('pymol OK')"))
if (Test-Path (Join-Path $EnvScripts "vina.exe")) {
    & (Join-Path $EnvScripts "vina.exe") --version
}
[void](Invoke-EnvPython @("-m", "dockflow_core.cli", "info"))

Write-Host ""
Write-Ok "Done. Activate and start:"
Write-Host "    conda activate $EnvName"
Write-Host "    dockflow info              # environment report"
Write-Host "    dockflow gui               # desktop application"
Write-Host "    dockflow run --config examples\configs\hiv1_protease_example.yaml"
Write-Host ""
Write-Host "or simply double-click  run_dockflow.bat  in the repository root."
