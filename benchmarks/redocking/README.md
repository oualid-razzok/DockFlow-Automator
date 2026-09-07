# Redocking benchmark suite

The core scientific-validation artifact of DockFlow (audit item 7): a set
of classic protein–ligand complexes re-docked with the standard pipeline,
reporting experimental-vs-predicted pose RMSD, success rate, docking
scores and runtimes.

**Status: the harness is complete and locally validated on 1HVR (best
pose RMSD 1.08 Å, success).  The full 24-complex table requires roughly
2–4 CPU-hours at exhaustiveness 8 and is intended to be executed on a
workstation/cluster; `results/` is seeded with the locally produced
rows and grows as the suite is run.**

## Methodology

For every complex in `complexes.csv` (24 entries, all verified against
RCSB to contain exactly one organic ligand after excluding waters, ions
and sugars; sources: PDBbind core set / Astex diverse set literature and
the HIV-1 protease series used by the example configs):

1. the PDB structure is downloaded from RCSB;
2. the receptor is prepared with the default `auto` engine
   (Gasteiger charges, polar hydrogens, waters/metals removed with
   warnings recorded in the manifest);
3. the search space is the co-crystallized ligand pocket plus 4 Å
   padding (recorded as `gridbox.source = pocket:<ligand>` in the
   manifest);
4. the ligand is prepared from the RCSB chemical-component definition
   (ideal coordinates re-embedded by RDKit/Meeko);
5. AutoDock Vina docks it with a fixed seed (2026) and configurable
   exhaustiveness (default 8);
6. every pose is scored by **symmetry-tolerant heavy-atom Kabsch RMSD**
   to the co-crystallized pose extracted from the PDB
   (`crystal_rmsd`; see `docs/interaction_criteria.md`);
7. **success** = best pose RMSD ≤ 2.0 Å, the conventional redocking
   criterion.

Every number is reproducible: each run directory contains the full
`manifest.json` (stage audit trail, preparation decisions) and
`environment.json` (scientific-stack fingerprint).

## Running

```bash
# quick subset (minutes)
python benchmarks/redocking/run_benchmark.py --pdb-ids 1HVR,1STP --exhaustiveness 8

# full suite (~2-4 CPU-hours at exhaustiveness 8; ~2x that at 16)
python benchmarks/redocking/run_benchmark.py --exhaustiveness 8

# highest-fidelity (matches the example config)
python benchmarks/redocking/run_benchmark.py --exhaustiveness 16
```

Outputs land in `results/`:

* `redocking_results.csv` — one row per complex (the core artifact);
* `summary.md` — the success-rate table;
* `<pdb_id>_benchmark/` — the full pipeline run directories.

## Expected results (literature context)

With AutoDock Vina 1.2 and rigid-receptor redocking, typical success
rates (top pose ≤ 2 Å) on curated sets of this kind are **50–75 %**;
small rigid ligands (biotin/streptavidin, benzamidine/trypsin) usually
recover sub-Å, while flexible peptidomimetic HIV-PR inhibitors are the
hard cases.  DockFlow's own committed rows are in `results/`; the
reproducibility workflow re-validates the 1HVR entry on every change.

## Cross-checking against a conventional workflow (item 8)

`../mgltools_comparison/` contains the harness and methodology for
comparing DockFlow against the legacy MGLTools
(`prepare_receptor4.py`/`prepare_ligand4.py`) → `vina` CLI workflow on
the same complexes, including the known Meeko-vs-MGLTools
parameterisation differences.
