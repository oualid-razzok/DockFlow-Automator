# DockFlow vs conventional MGLTools workflow (peer item 8)

Goal: demonstrate how DockFlow's automated preparation compares with the
conventional AutoDock Vina workflow (MGLTools `prepare_receptor4.py` /
`prepare_ligand4.py` -> Vina) on the redocking benchmark set, and
document the systematic differences.

## Status: RUN on the full 24-complex set (2026-09)

MGLTools is **Python-2-only legacy software** (final release 2011) and
cannot run in DockFlow's supported Python >= 3.10 environments; it is
therefore **not run in CI**.  The comparison was executed out-of-CI with
the **bioconda `mgltools 1.5.7` package**, which bundles its own
Python 2.7 interpreter (no system Python 2 required):

```bash
# one-time setup: extract the bioconda package and fix its placeholder
# paths (scripts/fix_mgltools_prefix.py in the DockFlow dev sandbox)
curl -sL -o mgltools.tar.bz2 \
  "https://api.anaconda.org/download/bioconda/mgltools/1.5.7/linux-64/mgltools-1.5.7-0.tar.bz2"
mkdir mgltools && tar xjf mgltools.tar.bz2 -C mgltools
python /path/to/scripts/fix_mgltools_prefix.py   # rewrites conda build
                                                 # placeholder prefixes
python benchmarks/mgltools_comparison/compare_workflow.py \
    --mgltools /path/to/mgltools
```

Results (see `results/summary.md` for the full per-complex table):

* 24 complexes attempted; **21 completed in the MGLTools arm**;
  **3 genuine MGLTools 1.5.7 crashes** (receptor prep on 1E66, ligand
  prep on 1MQ6 and 1SL3 - tracebacks preserved in `results/`); failures
  are recorded with reasons, never dropped.
* Pose recovery (best RMSD <= 2.0 A): **DockFlow 16/21 vs MGLTools
  14/21** on the MGL-completed subset.
* Mean RMSD delta (MGLTools - DockFlow): **+0.45 A** (DockFlow recovers
  the crystal pose slightly better on average); mean affinity delta
  -0.09 kcal/mol (MGLTools scores marginally stronger); receptor
  non-hydrogen atom counts agree exactly on every complex.

## Protocol

Both arms dock with the **same AutoDock Vina** (python bindings), the
**same search space** (the redocking benchmark's `gridbox.txt`, derived
from the co-crystal ligand), the **same seed (2026)** and
**exhaustiveness (8)**, `num_modes 9`.  The only difference is the
preparation engine:

1. **DockFlow arm** - the committed redocking benchmark results: the
   same cleaned receptor (`receptor_clean.pdb`) through the openbabel
   engine (Gasteiger, polar hydrogens), the same RCSB ideal ligand SDF
   through RDKit/Meeko.
2. **MGLTools arm** - `prepare_receptor4.py -A hydrogens -U nphs` on the
   same `receptor_clean.pdb`, `prepare_ligand4.py` on the same ideal SDF
   (converted to mol2 for ADT).

Because inputs, search space and docking engine are identical, every
delta in `results/comparison.csv` is attributable to preparation
chemistry (hydrogen model, charge assignment, torsion trees), which is
exactly what this benchmark measures.

## Known systematic differences (literature + measured)

* **Meeko vs prepare_ligand4**: both merge non-polar hydrogens onto
  their heavy atoms and use Gasteiger charges - but their torsion-tree
  detection differs on a small set of bonds, which can shift scores by
  ~0.1-0.5 kcal/mol and occasionally pose selection.
* **Macrocycles**: Meeko opens macrocyclic rings for flexible-macrocycle
  docking (G0/CG0 glue pseudo-atoms; 1HVR/XK2 in this set); MGLTools
  docks the closed ring as-is.
* **OpenBabel vs prepare_receptor4**: hydrogen placement rules differ
  slightly; measured receptor non-hydrogen atom counts agree on all
  completed complexes.

## Running / resuming

```bash
python benchmarks/mgltools_comparison/compare_workflow.py \
    --mgltools /path/to/mgltools --max-seconds 540   # resumable chunks
python benchmarks/mgltools_comparison/compare_workflow.py --report-only
```

Outputs in `results/`: `mgl_arm_results.csv` (raw MGLTools arm rows),
`comparison.csv` (both arms + deltas), `summary.md` (published table),
`comparison_stats.json` (machine-readable), `<pdb_id>/` per-complex
MGLTools artifacts (receptor/ligand PDBQTs, docked poses).
