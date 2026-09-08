"""Dependency license audit (audit item 38).

Generates ``licenses.json`` and ``THIRD_PARTY_LICENSES.md`` from the
installed distribution metadata plus a manual verification table for
packages whose metadata omits the license.  Re-run whenever dependencies
change; the output is committed.

Usage::

    python scripts/license_audit.py
"""

from __future__ import annotations

import json
import sys
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# distribution -> (extras, role) for everything DockFlow can pull in
PACKAGES: dict[str, tuple[str, str]] = {
    "requests": ("core", "HTTP downloads (RCSB/PubChem/ZINC/UniProt)"),
    "numpy": ("core", "numeric kernels (distances, Kabsch RMSD)"),
    "PyYAML": ("core", "pipeline YAML config parsing"),
    "meeko": ("prep", "ligand PDBQT preparation (prepare_ligand4 successor)"),
    "rdkit": ("prep", "chemistry: sanitisation, embedding, Gasteiger, SDF I/O"),
    "gemmi": ("prep", "mmCIF support (runtime dependency of meeko)"),
    "pandas": ("prep", "dataframe support (runtime dependency of meeko)"),
    "openbabel-wheel": ("obabel", "OpenBabel python bindings: receptor prep engine"),
    "vina": ("engine", "AutoDock Vina python bindings (docking engine)"),
    "PyQt6": ("gui", "desktop GUI toolkit"),
    "pyqt6-sip": ("gui", "PyQt6 support library"),
    "matplotlib": ("viz", "fallback renderer + ROC curves"),
    "dimorphite_dl": ("prep (optional)", "protonation-state enumeration (opt-in)"),
    "pytest": ("test (dev)", "test framework"),
    "pytest-cov": ("test (dev)", "coverage reporting"),
    "ruff": ("test (dev)", "linter/formatter"),
    "hatchling": ("build (dev)", "build backend"),
    "pybind11": ("build (dev)", "C++ accelerator bindings headers"),
    "scikit-build-core": ("build (dev)", "C++ accelerator build backend"),
}

# Licenses verified from the projects' repositories where the wheel
# metadata omits them (checked against upstream LICENSE files).
MANUAL_LICENSES: dict[str, tuple[str, str]] = {
    "openbabel-wheel": ("GPL-2.0-or-later",
                        "OpenBabel is GPL-2.0+; the wheel redistributes compiled OpenBabel code"),
    "PyQt6": ("GPL-3.0-only OR commercial (Riverbank)",
              "dual-licensed; GPL-3 terms apply to the free use"),
    "pyqt6-sip": ("BSD-2-Clause", "Riverbank SIP"),
    "meeko": ("LGPL-2.1", "forlilab/Meeko LICENSE"),
    "vina": ("Apache-2.0", "ccsb-scripps/AutoDock-Vina LICENSE"),
    "matplotlib": ("PSF-based (Matplotlib license)", "BSD-compatible, non-copyleft"),
    "dimorphite_dl": ("Apache-2.0", "durrantlab/dimorphite-dl"),
    "pytest": ("MIT", "pytest-dev/pytest LICENSE"),
    "ruff": ("MIT", "astral-sh/ruff LICENSE"),
    "hatchling": ("Apache-2.0", "ofek/hatch"),
    "pybind11": ("BSD-3-Clause", "pybind/pybind11 LICENSE"),
    "scikit-build-core": ("Apache-2.0", "scikit-build/scikit-build-core"),
}

# External tools DockFlow drives via subprocess (not python deps).
EXTERNAL_TOOLS: list[dict[str, str]] = [
    {"tool": "AutoDock Vina CLI (1.2.7)", "license": "Apache-2.0",
     "role": "docking engine (CLI backend)", "invoked": "subprocess; not bundled"},
    {"tool": "Smina", "license": "Apache-2.0 (fork of Vina)",
     "role": "docking engine option", "invoked": "subprocess; not bundled"},
    {"tool": "GNINA", "license": "Apache-2.0",
     "role": "CNN scoring/docking backend (item 32)",
     "invoked": "subprocess; not bundled"},
    {"tool": "obabel CLI", "license": "GPL-2.0-or-later",
     "role": "receptor prep fallback engine",
     "invoked": "subprocess; not bundled"},
    {"tool": "open-source PyMOL", "license": "BSD-like (PyMOL open-source license)",
     "role": "ray-traced rendering", "invoked": "subprocess; not bundled"},
    {"tool": "MGLTools (comparison only)", "license": "custom (MGLTools)",
     "role": "benchmark comparison harness, never a runtime dependency",
     "invoked": "external python2 environment, opt-in"},
    {"tool": "3Dmol.js (CDN)", "license": "BSD-3-Clause",
     "role": "interactive 3D viewer script inside generated interactive.html",
     "invoked": "loaded from cdn.jsdelivr.net at view time; never bundled"},
]


def entry_for(distribution: str, extra: str, role: str) -> dict:
    try:
        meta = metadata.metadata(distribution)
        version = metadata.version(distribution)
        license_field = (meta.get("License") or "").strip()
        classifiers = [c for c in meta.get_all("Classifier", [])
                       if c.startswith("License ::")]
    except metadata.PackageNotFoundError:
        return {
            "distribution": distribution, "version": "not installed",
            "license": MANUAL_LICENSES.get(distribution, ("unknown", ""))[0],
            "extra": extra, "role": role, "bundled_in_binary": False,
        }
    if distribution in MANUAL_LICENSES:
        license_name, note = MANUAL_LICENSES[distribution]
        source = "verified from upstream repository"
    elif classifiers:
        license_name = classifiers[0].replace(
            "License :: OSI Approved :: ", "").replace("License ::", "").strip()
        note = license_field[:80]
        source = "wheel classifier"
    else:
        license_name = license_field[:60] or "unknown"
        note = ""
        source = "wheel License field"
    return {
        "distribution": distribution,
        "version": version,
        "license": license_name,
        "license_note": note,
        "license_source": source,
        "extra": extra,
        "role": role,
        "bundled_in_binary": False,
    }


def main() -> int:
    entries = [entry_for(dist, extra, role)
               for dist, (extra, role) in PACKAGES.items()]
    payload = {
        "$note": ("Generated by scripts/license_audit.py from installed "
                  "distribution metadata + manually verified upstream "
                  "licenses. bundled_in_binary refers to redistributable "
                  "DockFlow artifacts (wheels/sdist ship no third-party code)."),
        "dockflow_license": "MIT",
        "python_dependencies": entries,
        "external_tools": EXTERNAL_TOOLS,
    }
    json_path = ROOT / "licenses.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# THIRD-PARTY LICENSES",
        "",
        "DockFlow-Automator is MIT-licensed and **ships no third-party code "
        "in its wheels/sdist**; the table below enumerates everything its "
        "extras can pull in at runtime (plus the external executables it "
        "drives via subprocess).  Machine-readable version: "
        "[`licenses.json`](licenses.json) (regenerate with "
        "`python scripts/license_audit.py`).",
        "",
        "## Redistribution notes (the ones that matter)",
        "",
        "- **OpenBabel is GPL-2.0-or-later** - including when used through "
        "its Python bindings (`openbabel-wheel`).  Distributing DockFlow "
        "together with the OpenBabel bindings in one aggregate makes the "
        "aggregate effectively GPL-2.0 terms territory for that component; "
        "the dependency-free and RDKit engine paths avoid it entirely, "
        "and the `[obabel]` extra is optional for exactly this reason.",
        "- **PyQt6 is GPL-3.0-only or commercial** (Riverbank dual "
        "license); the GUI is optional (`[gui]` extra) and the CLI/API "
        "paths do not use it.",
        "- **Meeko is LGPL-2.1** - dynamically imported, never modified "
        "or embedded; LGPL obligations are limited to the library itself.",
        "- Everything else (Vina Apache-2.0, RDKit BSD-3, NumPy BSD, "
        "matplotlib PSF-based, requests Apache-2.0, PyYAML MIT, pybind11 "
        "BSD-3) is permissive and unproblematic.",
        "",
        "## Python dependencies (by extra)",
        "",
        "| distribution | version | license | extra | role |",
        "|---|---|---|---|---|",
    ]
    for entry in entries:
        lines.append(
            f"| {entry['distribution']} | {entry['version']} "
            f"| {entry['license']} | {entry['extra']} | {entry['role']} |"
        )
    lines += [
        "",
        "## External tools (subprocess, never bundled)",
        "",
        "| tool | license | role |",
        "|---|---|---|",
    ]
    for tool in EXTERNAL_TOOLS:
        lines.append(f"| {tool['tool']} | {tool['license']} | {tool['role']} |")
    lines.append("")
    (ROOT / "THIRD_PARTY_LICENSES.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {json_path} and THIRD_PARTY_LICENSES.md "
          f"({len(entries)} python deps, {len(EXTERNAL_TOOLS)} external tools)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
