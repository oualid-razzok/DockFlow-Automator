# DockFlow vs conventional MGLTools workflow (peer item 8)

Goal: demonstrate that DockFlow's automation reproduces the
conventional AutoDock Vina workflow (MGLTools `prepare_receptor4.py` /
`prepare_ligand4.py` -> `vina` CLI) on the redocking benchmark set, and
document any systematic differences.

## Status and honest limitation

MGLTools is **Python-2-only legacy software** (final release 2011) and
cannot run in DockFlow's supported Python >= 3.10 environments; it is
therefore **not run in CI**.  The comparison protocol below is fully
specified and `compare_workflow.py` automates it whenever a MGLTools
environment is available (e.g. a conda env with `mgltools` from the
`binstar` channel, or the `mgltools` docker image).  The
preparation-level comparison (item 9) is performed with the modern
engines DockFlow actually ships and documented in
`docs/preparation_validation.md`.

## Protocol

For every complex in `benchmarks/redocking/complexes.csv`:

1. **Conventional arm**: `prepare_receptor4.py -r <pdb> -A hydrogens -U nphs`
   and `prepare_ligand4.py -l <ideal ligand>` (MGLTools, default
   Gasteiger charges and AD4 typing), then `vina --config` with the same
   box, seed and exhaustiveness as DockFlow.
2. **DockFlow arm**: the standard `dockflow run` pipeline (Meeko ligand
   preparation, OpenBabel/RDKit receptor preparation).
3. Compare per complex: best affinity (kcal/mol), best symmetry-aware
   RMSD, pose overlap (RMSD between the two workflows' best poses).

## Known systematic differences to expect (literature + item 9)

* **Meeko vs prepare_ligand4**: both merge non-polar hydrogens onto their
  heavy atoms and transfer the charge, and both use Gasteiger charges -
  but Meeko's torsion-tree detection differs in a small set of bonds
  (e.g. some amides), which can shift rotatable-bond counts and
  therefore scores by ~0.1-0.5 kcal/mol and occasionally pose selection.
* **OpenBabel/RDKit vs prepare_receptor4**: hydrogen placement rules and
  Gasteiger implementations differ slightly; atom typing for the
  receptor rarely differs because the AD4 type set is shared.
* **Ideal-SDF ligand geometry**: both arms use the same input, so the
  comparison isolates parameterisation, not embedding.

## Running (when MGLTools is available)

```bash
# with $MGLTOOLS/bin on PATH (python2 environment)
python benchmarks/mgltools_comparison/compare_workflow.py \
    --complexes ../redocking/complexes.csv --exhaustiveness 8 \
    --out comparison_results
```

Outputs `comparison.csv` (per-complex deltas) and a summary table for
`docs/preparation_validation.md`.
