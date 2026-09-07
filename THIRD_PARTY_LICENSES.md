# THIRD-PARTY LICENSES

DockFlow-Automator is MIT-licensed and **ships no third-party code in its wheels/sdist**; the table below enumerates everything its extras can pull in at runtime (plus the external executables it drives via subprocess).  Machine-readable version: [`licenses.json`](licenses.json) (regenerate with `python scripts/license_audit.py`).

## Redistribution notes (the ones that matter)

- **OpenBabel is GPL-2.0-or-later** - including when used through its Python bindings (`openbabel-wheel`).  Distributing DockFlow together with the OpenBabel bindings in one aggregate makes the aggregate effectively GPL-2.0 terms territory for that component; the dependency-free and RDKit engine paths avoid it entirely, and the `[obabel]` extra is optional for exactly this reason.
- **PyQt6 is GPL-3.0-only or commercial** (Riverbank dual license); the GUI is optional (`[gui]` extra) and the CLI/API paths do not use it.
- **Meeko is LGPL-2.1** - dynamically imported, never modified or embedded; LGPL obligations are limited to the library itself.
- Everything else (Vina Apache-2.0, RDKit BSD-3, NumPy BSD, matplotlib PSF-based, requests Apache-2.0, PyYAML MIT, pybind11 BSD-3) is permissive and unproblematic.

## Python dependencies (by extra)

| distribution | version | license | extra | role |
|---|---|---|---|---|
| requests | 2.32.5 | Apache Software License | core | HTTP downloads (RCSB/PubChem/ZINC/UniProt) |
| numpy | 2.1.3 | BSD License | core | numeric kernels (distances, Kabsch RMSD) |
| PyYAML | 6.0.3 | MIT License | core | pipeline YAML config parsing |
| meeko | 0.8.0 | LGPL-2.1 | prep | ligand PDBQT preparation (prepare_ligand4 successor) |
| rdkit | 2025.9.6 | BSD-3-Clause | prep | chemistry: sanitisation, embedding, Gasteiger, SDF I/O |
| gemmi | 0.7.5 | Mozilla Public License 2.0 (MPL 2.0) | prep | mmCIF support (runtime dependency of meeko) |
| pandas | 2.2.3 | BSD License | prep | dataframe support (runtime dependency of meeko) |
| openbabel-wheel | 3.1.1.23 | GPL-2.0-or-later | obabel | OpenBabel python bindings: receptor prep engine |
| vina | 1.2.7 | Apache-2.0 | engine | AutoDock Vina python bindings (docking engine) |
| PyQt6 | 6.11.0 | GPL-3.0-only OR commercial (Riverbank) | gui | desktop GUI toolkit |
| pyqt6-sip | 13.12.0 | BSD-2-Clause | gui | PyQt6 support library |
| matplotlib | 3.9.2 | PSF-based (Matplotlib license) | viz | fallback renderer + ROC curves |
| dimorphite_dl | 2.0.2 | Apache-2.0 | prep (optional) | protonation-state enumeration (opt-in) |
| pytest | 9.0.2 | MIT | test (dev) | test framework |
| pytest-cov | 7.0.0 | MIT License | test (dev) | coverage reporting |
| ruff | 0.16.6 | MIT | test (dev) | linter/formatter |
| hatchling | not installed | Apache-2.0 | build (dev) | build backend |
| pybind11 | not installed | BSD-3-Clause | build (dev) | C++ accelerator bindings headers |
| scikit-build-core | not installed | Apache-2.0 | build (dev) | C++ accelerator build backend |

## External tools (subprocess, never bundled)

| tool | license | role |
|---|---|---|
| AutoDock Vina CLI (1.2.7) | Apache-2.0 | docking engine (CLI backend) |
| Smina | Apache-2.0 (fork of Vina) | docking engine option |
| GNINA | Apache-2.0 | CNN scoring/docking backend (item 32) |
| obabel CLI | GPL-2.0-or-later | receptor prep fallback engine |
| open-source PyMOL | BSD-like (PyMOL open-source license) | ray-traced rendering |
| MGLTools (comparison only) | custom (MGLTools) | benchmark comparison harness, never a runtime dependency |
