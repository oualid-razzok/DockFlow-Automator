# Redocking benchmark suite

The core scientific-validation artifact of DockFlow (audit item 7): a set
of classic protein–ligand complexes re-docked with the standard pipeline,
reporting experimental-vs-predicted pose RMSD, success rate, docking
scores and runtimes.

**Status: RUN on all 24 complexes (2026-09, exhaustiveness 8, seed 2026,
cpu pinned for regeneration).  Results in `results/`: 24/24 completed;
pose recovery (FORMAL definition: rank-1 pose, symmetry-aware heavy-atom
RMSD ≤ 2.0 Å) 21/24 = 87.5 % (95 % Wilson CI 69.0–95.7 %); top-3
definition 23/24; legacy best-any-rank 23/24 = 95.8 % (CI 79.8–99.3 %);
mean best-pose RMSD 0.67 Å / median 0.62 Å; all definitions, thresholds,
sub-tables and failure-mode classification in `results/summary.md`.**

Protocol notes from the actual run:

* the published numbers use the atom-order-independent symmetry RMSD
  (corrected in v1.0.0-rc1 after the spyrmsd cross-check found the old
  implementation over-penalised reordered symmetric ligands; see
  `rmsd_crosscheck.py` and the CHANGELOG).  Under the corrected RMSD,
  only 1HTG fails at every rank (sampling failure: a 17-rotatable-bond
  peptidomimetic with no pose within 3.0 Å), while 1DWD and 1HPV fail
  the rank-1 definition only because their ≤ 2.0 Å pose ranks 3rd
  (ranking failures - exactly what the published top-N curve
  distinguishes);
* every failed/struggling complex keeps its full run directory
  (including all pose PDBQTs/SDFs) for inspection;
* complexes whose ligand appears in two copies (homodimers: 1HTG, 1Q41,
  1XLZ, 1N46) use ONE instance (the largest) for the box and the RMSD
  reference - the union of both copies would inflate the box past Vina's
  30 Å cap and break the atom-count match for the RMSD;
* every row is reproducible: fixed seed, deterministic Vina (thread
  count pinned by `regenerate.py`), recorded environment fingerprint;
  the full suite was regenerated and ASSERTED to match the committed
  CSV within 1e-6 (0 mismatches) - see `regenerate.py`.

## Regenerating / validating

```bash
python benchmarks/redocking/regenerate.py            # re-dock all 24 + assert match (item 14)
python benchmarks/redocking/regenerate.py --skip-docking   # dataset SHA-256 pins only
python benchmarks/redocking/run_benchmark.py --report-only  # summary/figures from the committed artifacts
python benchmarks/redocking/reanalyze.py             # re-run the ANALYSIS stage only
python benchmarks/redocking/rmsd_crosscheck.py       # spyrmsd cross-check (item 5)
```

## Methodology

For every complex in `complexes.csv` (24 entries, all verified against
RCSB to contain exactly one organic ligand after excluding waters, ions
and sugars; sources: PDBbind core set / Astex diverse set literature and
the HIV-1 protease series used by the example configs; the table pins
the ligand code, UniProt, resolution and the SHA-256 of the committed
raw PDB):

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
6. every pose is scored by **symmetry-aware heavy-atom RMSD**
   (graph-automorphism minimised, atom-order independent) to the
   co-crystallized pose extracted from the PDB
   (`crystal_rmsd`; see `docs/interaction_criteria.md`);
7. **success** is published under EVERY definition (rank-1 / top-3 /
   best-any-rank at 1.0/2.0/3.0 Å thresholds) so the number cannot be
   definition-shopped; the formal definition is rank-1 ≤ 2.0 Å.

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
