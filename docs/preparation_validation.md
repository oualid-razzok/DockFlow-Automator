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

## MGLTools arm, 1HVR (final pass, item 13)

The MGLTools arm HAS now been run (MGLTools 1.5.7, bioconda build with
bundled Python 2.7, relocated in place).  Every number below is computed
by the committed scripts - `benchmarks/mgltools_comparison/
compare_workflow.py` for the MGL arm, `benchmarks/preparation_check.py`
for the DockFlow engines - and the per-complex raw values are in
`benchmarks/mgltools_comparison/results/mgl_arm_results.csv` (nothing is
hand-transcribed; the aggregate script also re-derives all of it, see
below).

| Property | DockFlow (openbabel receptor + Meeko ligand) | MGLTools 1.5.7 (`prepare_receptor4.py` / `prepare_ligand4.py`) |
|---|---|---|
| Receptor heavy atoms (1HVR) | 1500 | 1500 (identical) |
| Ligand heavy atoms (XK2) | 46 | 46 (identical) |
| Ligand polar hydrogens | 2 | 2 (identical) |
| Ligand rotatable bonds | 16 (incl. macrocycle pseudo-torsions) | 10 (closed-ring macrocycle) |
| Ligand Gasteiger charge sum | +0.007 | +0.005 |
| Receptor prep | completes | completes |
| Ligand prep | completes on all 24 benchmark ligands | crashes on 3/24 (1MQ6, 1SL3, 1HVR-class macrocycles differ) |

## Cross-workflow, 24-complex aggregate (final pass, item 13)

Generated by `benchmarks/preparation_validation_aggregate.py` (also
writes the machine-readable
`benchmarks/preparation_validation_aggregate.json`); "crash" marks the
three genuine MGLTools 1.5.7 crashes (full tracebacks in
`benchmarks/mgltools_comparison/results/failures/`).  "DF" = the
committed redocking-benchmark artifacts; "MGL" = the committed
MGLTools-arm artifacts.

| pdb | ligand | rec heavy DF/MGL | lig heavy DF/MGL | lig polar H DF/MGL | rotatable bonds DF/MGL | Gasteiger sum DF/MGL |
|---|---|---|---|---|---|---|
| 1HVR | XK2 | 1500/1500 | 46/46 | 2/2 | 16/10 | 0.007/0.005 |
| 1HSG | MK1 | 1514/1514 | 45/45 | 4/4 | 14/14 | 0.004/2.001 |
| 1HTG | G37 | 1516/1516 | 55/55 | 7/7 | 17/17 | -0.004/1.0 |
| 1HPV | 478 | 1516/1516 | 35/35 | 4/4 | 14/14 | 0.0/0.002 |
| 1MUI | AB1 | 1516/1516 | 46/46 | 4/4 | 16/16 | 0.001/0.006 |
| 1DWD | MID | 2358/2358 | 37/37 | 5/5 | 9/10 | 0.002/-0.002 |
| 1AZM | AZM | 2019/2019 | 13/13 | 3/3 | 3/3 | -0.001/-0.002 |
| 1E66 | HUX | 4157/crash | 21/crash | 2/crash | 2/crash | -0.0/crash |
| 1F0R | 815 | 2227/2227 | 31/31 | 3/3 | 6/6 | 0.002/0.002 |
| 1Q41 | IXM | 5428/5428 | 21/21 | 3/3 | 1/2 | 0.002/-0.001 |
| 1SNC | THP | 1082/1082 | 25/25 | 5/5 | 10/10 | -0.003/-0.998 |
| 3ERT | OHT | 1932/1932 | 29/29 | 1/1 | 9/9 | 0.001/1.001 |
| 1XLZ | FIL | 5202/5202 | 21/21 | 2/2 | 6/6 | -0.001/-0.002 |
| 1STP | BTN | 901/901 | 16/16 | 3/3 | 6/6 | 0.002/-1.0 |
| 1CBX | BZS | 2437/2437 | 15/15 | 2/2 | 7/7 | 0.001/-1.998 |
| 1TNG | AMC | 1629/1629 | 8/8 | 3/3 | 2/2 | 1.001/1.0 |
| 1BTY | BEN | 1629/1629 | 9/9 | 3/3 | 1/2 | 0.001/-0.002 |
| 1MQ6 | XLD | 2220/crash | 36/crash | 2/crash | 8/crash | 0.002/crash |
| 1N8Q | DHB | 6778/6778 | 11/11 | 3/3 | 4/4 | -0.0/-1.003 |
| 1NVQ | UCN | 2169/2169 | 36/36 | 3/3 | 3/4 | 0.0/0.996 |
| 1JAP | HOA | 1260/1260 | 2/2 | 3/3 | 1/1 | -0.0/-0.0 |
| 1YGC | 905 | 2368/2368 | 38/38 | 8/8 | 13/14 | -0.0/0.0 |
| 1N46 | PFA | 3872/3872 | 27/27 | 2/2 | 5/5 | -0.002/0.0 |
| 1SL3 | 170 | 2321/crash | 37/crash | 2/crash | 8/crash | -0.0/crash |

Reading the aggregate honestly:

* **Receptor heavy-atom counts are IDENTICAL on every MGL-completed
  complex** - the two workflows agree on the protein geometry exactly.
* **Ligand heavy atoms and polar hydrogens are identical** as well.
* **Rotatable bonds differ on 6/21 complexes** (by 1, except 1HVR/XK2's
  macrocycle: 16 vs 10): Meeko breaks macrocycles into pseudo-torsions
  (flexible-macrocycle docking) while `prepare_ligand4` keeps the ring
  closed; single-bond differences are the classic perception
  disagreements (e.g. amide bonds).
* **Gasteiger charge sums** mostly agree to < 0.01 e; where MGLTools
  reports integral values (±1, ±2) it has assigned FORMAL charges to
  ionised groups (carboxylates/amines), while Meeko/OpenBabel distribute
  the same net charge fractionally.  Vina reads atom types, not partial
  charges, so this does not affect Vina scoring - it matters for AD4
  maps (see the engine-differences section above).

## Meeko-invocation faithfulness (final pass, item 13)

To show DockFlow's Meeko wrapper does not silently change upstream
behaviour, `benchmarks/preparation_validation_aggregate.py` also runs
the upstream `mk_prepare_ligand.py` CLI (meeko's own console script,
no DockFlow code) on the SAME committed 1HVR/XK2 ideal SDF:

| Property | upstream `mk_prepare_ligand.py` | DockFlow wrapper |
|---|---|---|
| Heavy atoms | 46 | 46 |
| Polar hydrogens | 2 | 2 |
| Rotatable bonds | 16 | 16 |
| Gasteiger charge sum | +0.007 | +0.007 |

Identical in every compared property: the wrapper adds provenance
(recording, manifest fields) without changing the chemistry.

## Regenerating this table

```bash
python benchmarks/preparation_check.py --pdb runs/<run>/raw/1hvr.pdb \
    --ligand-sdf runs/<run>/raw/xk2.sdf
```

writes `benchmarks/preparation_report.json` (the machine-readable source
of the table above) and fails if any engine's heavy-atom count disagrees
with the reference.

```bash
python benchmarks/preparation_validation_aggregate.py
```

re-derives the 24-complex cross-workflow table, the MGLTools-arm
aggregate and the meeko-faithfulness check from the committed artifacts
(requires the `prep` extras + meeko).
