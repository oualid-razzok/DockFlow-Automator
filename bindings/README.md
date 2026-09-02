# dockflow-bindings

Optional C++ accelerator kernels for **DockFlow-Automator**, built with
pybind11 + scikit-build-core.  The Python package works without them (pure
NumPy fallbacks produce identical results); install them for faster batch
analysis of large pose libraries.

## Functions

| function | purpose |
|---|---|
| `parse_pdbqt_atoms(text)` | fast PDBQT `ATOM`/`HETATM` parsing |
| `grid_box(coords, padding, min_size)` | bounding box + padding |
| `pairwise_min_dist(a, b)` | per-atom minimum distance between clouds |
| `min_contacts(a, b, cutoff)` | all pairs within a cutoff (contact map) |
| `direct_rmsd(a, b)` | RMSD without superposition |
| `kabsch_rmsd(a, b)` | RMSD after optimal superposition (Horn quaternion method, in-house Jacobi eigensolver - no BLAS dependency) |
| `ligand_efficiency(affinity, heavy)` | affinity per heavy atom |
| `box_corners(center, size)` | 8 grid-box corners |

## Build

```bash
# from the repository root
pip install ./bindings

# or with a plain CMake workflow
cmake -B bindings/build -S bindings -DPython3_EXECUTABLE=$(which python)
cmake --build bindings/build -j
cmake --install bindings/build
```

Requires a C++17 compiler and CMake >= 3.18.  pybind11 is fetched
automatically by CMake when not installed.
