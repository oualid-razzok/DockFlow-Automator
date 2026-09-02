"""Shared pytest fixtures: embedded structures, ligands, logs, CLI hooks."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make the repository importable without installation.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# ---------------------------------------------------------------------------
# A small but realistic receptor: 6-residue peptide + ZN + water + ligand
# ---------------------------------------------------------------------------
RECEPTOR_PDB = """\
HEADER    TEST    DockFlow-Automator fixtures
COMPND    MINI PROTEIN FOR DOCKFLOW TESTS
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.006   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.216   2.403   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       2.040  -0.770  -1.200  1.00 20.00           C
ATOM      6  N   GLY A   2       3.476   1.420   0.000  1.00 20.00           N
ATOM      7  CA  GLY A   2       4.024   2.840   0.000  1.00 20.00           C
ATOM      8  C   GLY A   2       5.494   2.840   0.000  1.00 20.00           C
ATOM      9  O   GLY A   2       6.284   1.857   0.000  1.00 20.00           O
ATOM     10  N   SER A   3       6.042   4.260   0.000  1.00 20.00           N
ATOM     11  CA  SER A   3       7.511   4.260   0.000  1.00 20.00           C
ATOM     12  C   SER A   3       8.059   5.680   0.000  1.00 20.00           C
ATOM     13  O   SER A   3       7.269   6.663   0.000  1.00 20.00           O
ATOM     14  CB  SER A   3       8.093   3.490   1.200  1.00 20.00           C
ATOM     15  OG  SER A   3       9.479   3.490   1.200  1.00 20.00           O
ATOM     16  N   LYS A   4       9.529   5.680   0.000  1.00 20.00           N
ATOM     17  CA  LYS A   4      10.078   7.100   0.000  1.00 20.00           C
ATOM     18  C   LYS A   4      11.548   7.100   0.000  1.00 20.00           C
ATOM     19  O   LYS A   4      12.338   6.117   0.000  1.00 20.00           O
ATOM     20  CB  LYS A   4      10.446   7.870  -1.200  1.00 20.00           C
ATOM     21  CG  LYS A   4      11.916   7.870  -1.200  1.00 20.00           C
ATOM     22  CD  LYS A   4      12.386   8.000   0.200  1.00 20.00           C
ATOM     23  CE  LYS A   4      13.856   8.000   0.200  1.00 20.00           C
ATOM     24  NZ  LYS A   4      14.326   8.130   1.600  1.00 20.00           N
ATOM     25  N   GLU A   5      12.096   8.520   0.000  1.00 20.00           N
ATOM     26  CA  GLU A   5      13.566   8.520   0.000  1.00 20.00           C
ATOM     27  C   GLU A   5      14.114   9.940   0.000  1.00 20.00           C
ATOM     28  O   GLU A   5      13.324  10.923   0.000  1.00 20.00           O
ATOM     29  CB  GLU A   5      14.196   7.750  -1.200  1.00 20.00           C
ATOM     30  CG  GLU A   5      15.666   7.750  -1.200  1.00 20.00           C
ATOM     31  CD  GLU A   5      16.136   7.880   0.200  1.00 20.00           C
ATOM     32  OE1 GLU A   5      17.426   7.880   0.200  1.00 20.00           O
ATOM     33  OE2 GLU A   5      15.516   8.010   1.500  1.00 20.00           O
ATOM     34  N   VAL A   6      15.584   9.940   0.000  1.00 20.00           N
ATOM     35  CA  VAL A   6      16.132  11.360   0.000  1.00 20.00           C
ATOM     36  C   VAL A   6      17.602  11.360   0.000  1.00 20.00           C
ATOM     37  O   VAL A   6      18.392  10.377   0.000  1.00 20.00           O
ATOM     38  CB  VAL A   6      15.762  12.130   1.200  1.00 20.00           C
ATOM     39  CG1 VAL A   6      16.232  12.260   2.600  1.00 20.00           C
ATOM     40  CG2 VAL A   6      14.292  12.130   1.200  1.00 20.00           C
HETATM   41  ZN  ZN2 A  99       9.000   9.000   9.000  1.00 30.00          ZN
HETATM   42  O   HOH A 100      10.000  10.000  10.000  1.00 30.00           O
HETATM   43  O   HOH A 101      11.000  10.000  10.000  0.50 30.00           O
HETATM   44  C1  BEN A 901      12.000   9.000   8.000  1.00 25.00           C
HETATM   45  C2  BEN A 901      13.300   9.500   8.400  1.00 25.00           C
HETATM   46  C3  BEN A 901      13.700  10.800   8.000  1.00 25.00           C
HETATM   47  C4  BEN A 901      12.900  11.700   7.300  1.00 25.00           C
HETATM   48  C5  BEN A 901      11.600  11.200   6.900  1.00 25.00           C
HETATM   49  C6  BEN A 901      11.200   9.900   7.300  1.00 25.00           C
HETATM   50  O1  BEN A 901      10.500   9.400   6.500  0.50 25.00           O
HETATM   50A O1  BEN A 901      10.600   9.500   6.400  1.00 25.00           O
TER
END
"""

# A minimal valid PDBQT ligand (methanol-ish) with 2 torsions.
LIGAND_PDBQT = """\
REMARK  4 active site residues for ligand LIG
REMARK  6 long-range interactions
ROOT
ATOM      1  C1  LIG A   1      12.500   9.100   8.100  1.00  0.00     0.214 C
ATOM      2  O1  LIG A   1      13.100   9.800   8.700  1.00  0.00    -0.641 OA
ATOM      3  H1  LIG A   1      12.900  10.700   9.100  1.00  0.00     0.427 HD
ENDROOT
BRANCH   1   2
ATOM      4  C2  LIG A   1      11.100   9.300   8.300  1.00  0.00     0.114 C
ATOM      5  H2  LIG A   1      10.600   8.400   8.000  1.00  0.00     0.031 H
ENDROOT
TORSDOF 2
"""

# A docked multi-pose output in Vina 1.2 style.
DOCKED_PDBQT = """\
MODEL 1
REMARK VINA RESULT:    -9.423      0.000      0.000
REMARK SMILES bogus
ROOT
ATOM      1  C1  LIG A   1      12.500   9.100   8.100  1.00  0.00     0.214 C
ATOM      2  O1  LIG A   1      13.100   9.800   8.700  1.00  0.00    -0.641 OA
ATOM      3  H1  LIG A   1      12.900  10.700   9.100  1.00  0.00     0.427 HD
ENDROOT
TORSDOF 2
ENDMDL
MODEL 2
REMARK VINA RESULT:    -8.711      1.234      2.100
ROOT
ATOM      4  C1  LIG A   1      12.100   9.500   8.400  1.00  0.00     0.214 C
ATOM      5  O1  LIG A   1      12.700  10.200   9.000  1.00  0.00    -0.641 OA
ATOM      6  H1  LIG A   1      12.500  11.100   9.400  1.00  0.00     0.427 HD
ENDROOT
TORSDOF 2
ENDMDL
MODEL 3
REMARK VINA RESULT:    -7.905      3.412      4.823
ROOT
ATOM      7  C1  LIG A   1      13.500   8.500   7.800  1.00  0.00     0.214 C
ATOM      8  O1  LIG A   1      14.100   9.200   8.400  1.00  0.00    -0.641 OA
ATOM      9  H1  LIG A   1      13.900  10.100   8.800  1.00  0.00     0.427 HD
ENDROOT
TORSDOF 2
ENDMDL
"""

VINA_LOG = """\
AutoDock Vina v1.2.5
 #################################################################
 # If you used AutoDock Vina in your work, please cite:          #
 ...
 #################################################################

Computing Vina grid ... done.
Performing docking (random seed = 2026) ...

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b. | rmsd u.b.
-----+------------+----------+----------
   1        -9.423          0.000      0.000
   2        -8.711          1.234      2.100
   3        -7.905          3.412      4.823
"""

# A tiny valid SDF (methanol: C, O and 4 hydrogens).
METHANOL_SDF = """\
methanol
  DockFlow test fixture

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4500    0.9000    0.3000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4500   -0.9000    0.3000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4500    0.0000   -1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8000   -0.8000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
M  END
$$$$
"""


@pytest.fixture
def receptor_pdb_path(tmp_path: Path) -> Path:
    path = tmp_path / "mini_receptor.pdb"
    path.write_text(RECEPTOR_PDB, encoding="utf-8")
    return path


@pytest.fixture
def receptor_pdb_text() -> str:
    return RECEPTOR_PDB


@pytest.fixture
def ligand_pdbqt_path(tmp_path: Path) -> Path:
    path = tmp_path / "ligand.pdbqt"
    path.write_text(LIGAND_PDBQT, encoding="utf-8")
    return path


@pytest.fixture
def ligand_pdbqt_text() -> str:
    return LIGAND_PDBQT


@pytest.fixture
def docked_pdbqt_path(tmp_path: Path) -> Path:
    path = tmp_path / "ligand_out.pdbqt"
    path.write_text(DOCKED_PDBQT, encoding="utf-8")
    return path


@pytest.fixture
def docked_pdbqt_text() -> str:
    return DOCKED_PDBQT


@pytest.fixture
def vina_log_text() -> str:
    return VINA_LOG


@pytest.fixture
def methanol_sdf_path(tmp_path: Path) -> Path:
    path = tmp_path / "methanol.sdf"
    path.write_text(METHANOL_SDF, encoding="utf-8")
    return path


@pytest.fixture
def run_dir(tmp_path: Path) -> Path:
    directory = tmp_path / "run"
    directory.mkdir()
    return directory
