#!/usr/bin/env bash
# =============================================================================
# End-to-end example: HIV-1 protease redocking (PDB 1HVR / ligand XK2).
#
# Requires network access and at least the [prep,engine] extras:
#   pip install "dockflow-automator[prep,engine]"
#
# Usage:
#   bash examples/run_example.sh [workdir]
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKDIR="${1:-runs}"
CLI="python ${ROOT}/scripts/dockflow_cli.py"

say() { printf '\n\033[1;36m=== %s ===\033[0m\n' "$*"; }

say "1/6 downloading the target (PDB 1HVR)"
$CLI download pdb --id 1HVR --out "${WORKDIR}/raw"

say "2/6 downloading the co-crystallized ligand (XK2)"
$CLI download ligand --pdb-ligand XK2 --out "${WORKDIR}/raw"

say "3/6 preparing the receptor"
$CLI prep receptor --in "${WORKDIR}/raw/1hvr.pdb" --out-dir "${WORKDIR}/prepared"

say "4/6 preparing the ligand from RCSB ideal coordinates"
$CLI prep ligand --in "${WORKDIR}/raw/xk2.sdf" --name xk2 --out-dir "${WORKDIR}/prepared"

say "5/6 computing the grid box around the reference ligand"
$CLI gridbox --structure "${WORKDIR}/raw/1hvr.pdb" --resname XK2 \
    --padding 4 --out "${WORKDIR}/gridbox.txt"

say "6/6 docking (AutoDock Vina)"
$CLI dock \
    --receptor "${WORKDIR}/prepared/receptor.pdbqt" \
    --ligands "${WORKDIR}/prepared/xk2.pdbqt" \
    --config "${WORKDIR}/gridbox.txt" \
    --exhaustiveness 8 --num-modes 9 \
    --out-dir "${WORKDIR}/docking"

say "analysis"
$CLI analyze --docking "${WORKDIR}/docking" \
    --receptor "${WORKDIR}/prepared/receptor.pdbqt" --out-dir "${WORKDIR}/analysis"

say "done - inspect ${WORKDIR}/docking/summary.csv and ${WORKDIR}/analysis/"
