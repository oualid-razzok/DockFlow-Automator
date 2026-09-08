# Enrichment benchmark: actives vs decoys (audit item 22; final pass items 10-12)

Quantifies the "decoys score worse" claim with a standardised
active/decoy panel instead of arbitrary aspirin/caffeine controls.
Every number in `results/` is regenerable from the committed scripts +
committed datasets (see "Reproduce" at the end).

## Status: RUN on the real DUD-E panel - 3 targets (final pass, item 12)

| target | receptor (PDB) | n actives | n decoys | ROC AUC | EF@1% | EF@5% | BEDROC (α=20) |
|---|---|---|---|---|---|---|---|
| hivpr | 1HVR | 40 | 93 | 0.665 | 3.325 | 2.375 | 0.6647 |
| aa2ar | 3EML | 40 | 93 | 0.593 | 1.663 | 1.900 | 0.5486 |
| parp1 | 2RD6 | 40 | 93 | 0.729 | 3.325 | 1.425 | 0.5157 |

Per-target summaries: `results/<target>/summary.md`; multi-target table:
`results/summary.md`.  This is an **evaluation** of the DockFlow+Vina
pipeline, NOT a claim of superior virtual-screening performance
(single-seed, exhaustiveness 1, no comparative baseline, 3 of 102 DUD-E
targets - see each summary's caveats).

## Panel: DUD-E (committed, SHA-256 pinned)

The DUD-E ISM lists are COMMITTED under `data/` so the exact panels are
reproducible without re-downloading (the SHA-256 pins let a reviewer
detect upstream drift):

| file | SHA-256 (first 16 hex) |
|---|---|
| `data/hivpr_actives_final.ism` | `60f6b3ccbb815f10` (536 actives) |
| `data/hivpr_decoys_final.ism` | `bf67ae7fcec2f822` (35 750 decoys) |
| `data/aa2ar_actives_final.ism` | `690845896b1f1c0b` (482 actives) |
| `data/aa2ar_decoys_final.ism` | `b572781d1316c3b1` (31 550 decoys) |
| `data/parp1_actives_final.ism` | `78eeb68d302a118e` (508 actives) |
| `data/parp1_decoys_final.ism` | `39b177e9e3016cb3` (30 050 decoys) |

Source URLs: `https://dude.docking.org/targets/<target>/actives_final.ism`
and `.../decoys_final.ism` (no registration).  Full digests are
recomputable with `sha256sum data/*.ism`.

## Subsample (exact, seeded)

Each panel is a seeded subsample of 40 actives + 93 decoys:

* RNG: `random.Random(seed).sample(list, n)` (Python stdlib,
  deterministic across platforms for a fixed seed);
* seed: **2026** (the same seed as the docking);
* the ISM files are read line-by-line in file order (pinned by the
  SHA-256 above), so `sample` sees the same input list everywhere;
* verified: the committed `results/hivpr/enrichment_results.csv`
  ligand-id set is byte-for-byte the reproduction of
  `random.Random(2026).sample` over the committed ISM files (133/133
  identical ids).

## Preparation (ligand + receptor)

* **Ligand prep**: Meeko via `LigandPreparator(LigandPrepOptions(
  embed_3d=True, minimize=True, random_seed=2026))` — RDKit 3D
  embedding + MMFF94 minimisation (seeded, deterministic), Gasteiger
  charges, non-polar hydrogens merged, AD4 types, macrocycle breaking.
* **Protonation state**: INPUT state, no pH adjustment (no
  dimorphite-dl) — recorded per ligand in the run manifest.
* **Tautomers**: none enumerated (single input tautomer).
* **Salt stripping**: Meeko default (salts and counterions removed).
* **Receptor prep**: DockFlow `ReceptorPreparator` (engine `auto` →
  openbabel in the published run), Gasteiger charges, deposited
  protonation + template hydrogens, waters/metals/cofactors removed
  with per-species counts recorded.
* **hivpr protocol note**: the receptor PDBQT and grid box are reused
  VERBATIM from the redocking benchmark's committed 1HVR run
  (`benchmarks/redocking/results/1HVR_benchmark/`) so the enrichment
  run is bit-identical to the published one; aa2ar (3EML) and parp1
  (2RD6) are downloaded from RCSB and prepared by the same code path.

## Grid box (exact)

Ligand-derived (co-crystal ligand of the reference structure, 4.0 Å
padding, assumption strength: strong):

| target | reference PDB | box ligand | gridbox.txt |
|---|---|---|---|
| hivpr | 1HVR | XK2 | `results/hivpr/gridbox.txt` (copied from the redocking benchmark) |
| aa2ar | 3EML | ZMA | `results/aa2ar/gridbox.txt` (generated: ligand + 4.0 Å) |
| parp1 | 2RD6 | 78P | `results/parp1/gridbox.txt` (generated: ligand + 4.0 Å) |

## Docking (exact)

* Backend: AutoDock Vina **python bindings** (VinaConfig, backend
  `"python"`), Vina **1.2.7**;
* scoring: `vina`;
* exhaustiveness: **1**, num_modes: **1**, seed: **2026**, cpu: 0
  (all available threads - see the reproducibility caveat below);
* per-ligand timeout: 600 s.

## Scoring and tie-breaking

* Ranking score: the single Vina output **affinity** column
  (`best_affinity_kcal_mol`; lower = better; negated internally so
  higher = better ranked).
* Ties in score are broken by **ascending ligand identifier** (the
  ligand list is sorted before metric computation) - deterministic.
* Duplicate handling: DUD-E distributes already-deduplicated lists;
  this is VERIFIED, not trusted - the runner reports duplicate
  canonical-SMILES entries per panel (hivpr: 7 duplicate entries,
  reported in `enrichment_metrics.json` → `duplicates`).  Duplicates
  are never dropped: the published panel is exactly the seeded
  subsample.

## Metrics (exact formulas)

All implemented in `dockflow_core/enrichment.py` (verified against
RDKit's Truchon-Bayly reference implementation to 1e-12):

* **ROC AUC**: rank-based (Mann-Whitney U / n_actives · n_decoys).
  Random baseline: 0.5.
* **EF@x%**: (actives in top x% of the ranking / (x% of N)) / actives
  fraction — enrichment factor at x% = 1% and 5%.  Random baseline: 1.0.
* **BEDROC** (Truchon & Bayly 2007, α = 20):
  `RIE = Σ_actives exp(-α·rank_i/N) / (n_act · denom)`,
  `denom = (1 - exp(-α)) / (N · (exp(α/N) - 1))`, then
  BEDROC = (RIE − RIEmin)/(RIEmax − RIEmin) with
  `RIEmax = (1 − exp(-α·R)) / (R·(1 − exp(-α)))`,
  `RIEmin = (1 − exp(α·R)) / (R·(1 − exp(α)))`, R = n_actives/N.
  The random baseline is NOT 0.5 at α = 20 — it depends on the actives
  fraction, so the runner MEASURES it: mean over 10 000 seeded random
  permutations of the same panel (recorded per target in
  `enrichment_metrics.json` → `random_baselines`).

## Reproduce

Re-run the exact published hivpr experiment and assert the metrics match
within 1e-6:

```bash
python benchmarks/enrichment/run_enrichment.py --reproduce
```

Reproduces: the committed ISM panel → the seeded 133-ligand subsample →
Meeko prep (same seed) → Vina docking (same seed) → metric comparison
against `results/hivpr/enrichment_metrics.json`.  Exit code 0 = match.

**Thread caveat (honest):** Vina's search is thread-count dependent
(the v0.3.1 CI lesson).  The published runs used cpu=0 on a 2-vCPU
machine; on the same thread count the metrics reproduce exactly; on a
different thread count the search trajectories differ and `--reproduce`
will report the mismatch instead of silently passing.  Pin the thread
count (run with the same `cpu`) when reproducing on other hardware.

Other entry points:

```bash
python benchmarks/enrichment/run_enrichment.py --target hivpr   # resumable
python benchmarks/enrichment/run_enrichment.py --target aa2ar
python benchmarks/enrichment/run_enrichment.py --target parp1
python benchmarks/enrichment/run_enrichment.py --report-only --target hivpr
```
