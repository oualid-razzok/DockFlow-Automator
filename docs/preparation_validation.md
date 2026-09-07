# Preparation validation: cross-engine comparison (peer item 9)

This document compares DockFlow's preparation outputs across its
interchangeable engines on the example target (HIV-1 protease, PDB
**1HVR**, ligand **XK2**), so users know exactly which scientific
decisions differ between engines.  All numbers below are **computed by
the committed script** `benchmarks/preparation_check.py` and regenerated
on demand; nothing is transcribed by hand.

## Receptor preparation (PDB 1HVR, default options)

| Property | OpenBabel engine | RDKit engine | `none` engine (fallback) |
|---|---|---|---|
| Total atoms out | 1846 | 1848 | 1826 |
| Heavy atoms | 1500 | 1500 | 1500 |
| Polar hydrogens (HD) | 346 | 348 | 326 (from the PDB when present) |
| Gasteiger charge sum | +0.07 | **+23.97** | 0 (no charges) |
| Charge range | −0.51 … +0.39 | −0.51 … +0.40 | 0 |
| AD4 types: aromatic C (A) | 74 | 74 | 0 (no aromaticity perception) |
| AD4 types: N / NA | 256 / 2 | 248 / 10 | 258 / 0 |
| AD4 types: OA | 268 | 268 | 268 |
| HETATM removed (recorded) | CSO 18, XK2 46 | same | same |
| Warnings | 1 (removed species) | 1 | 2 |

### Interpretation (the scientifically relevant differences)

1. **Heavy-atom counts agree exactly** (1500) across all engines — the
   protein geometry is identical; only hydrogen treatment differs.
2. **Charge model**: OpenBabel's Gasteiger assignment distributes
   residual partial charges that sum to ~0, while RDKit keeps the
   residue-template **formal charges** (net +23.97 for this structure's
   protonation assignment).  *For Vina and Vinardo scoring this is
   irrelevant* — those scoring functions read AutoDock atom **types**
   and geometry, not partial charges.  For `scoring: ad4` (precomputed
   maps) partial charges matter directly; use one engine consistently
   and record it (the manifest does).
3. **Aromatic typing**: the `none` engine cannot perceive aromaticity
   (all carbons typed `C`, no `A`), which affects hydrophobic contact
   classification and AD4 maps, but not Vina geometry scoring.
4. **NA/OA acceptor typing** differs slightly (2 vs 10 `NA`) because the
   engines disagree about amide nitrogen tautomer/acceptor status for a
   few residues — bounded by the shared typing rules after the graph is
   perceived.
5. The `none` engine adds no hydrogens and assigns zero charges: its
   output is a **degraded geometry-only receptor** (see the per-engine
   "Known scientific differences" section in the README), suitable for
   CI testing, not for science.

## Ligand preparation (XK2 from the RCSB ideal SDF, Meeko)

| Property | Value |
|---|---|
| Real heavy atoms | 46 (41 C, 3 O, 2 N — matches the crystal copy exactly) |
| Meeko virtual atoms (G0/CG0, guanidinium typing) | 2 (excluded from RMSD computations) |
| Polar hydrogens (HD) | 2 |
| Non-polar hydrogens | merged into heavy atoms (charge transferred) |
| Rotatable bonds | 10 |
| Charge sum | +0.007 |
| Input hydrogens | none in the RCSB ideal SDF → all added by RDKit/Meeko (recorded per ligand in `manifest.ligands[*].prep`) |

Meeko is the modern successor of `prepare_ligand4.py` and reproduces its
core behaviour (Gasteiger charges, non-polar-H merging, torsion tree);
known differences are documented in
`benchmarks/mgltools_comparison/README.md`.

## MGLTools comparison protocol (peer items 8/9)

The legacy MGLTools scripts are Python-2-only and cannot run in
DockFlow's Python ≥ 3.10 environments; the automated comparison harness
lives in `benchmarks/mgltools_comparison/` and is executed whenever a
MGLTools environment is available.  The protocol, the comparison columns
(atom counts, partial-charge statistics, AutoDock type histograms,
rotatable-bond counts, protonation states) and the expected systematic
differences are fully specified there; this table will be extended with
the MGLTools arm's numbers when that environment is used.

## Regenerating this table

```bash
python benchmarks/preparation_check.py --pdb runs/<run>/raw/1hvr.pdb \
    --ligand-sdf runs/<run>/raw/xk2.sdf
```

writes `benchmarks/preparation_report.json` (the machine-readable source
of the table above) and fails if any engine's heavy-atom count disagrees
with the reference.
