# DUD-E enrichment evaluation (target hivpr, subsample 133 ligands)

Enrichment evaluation on the DUD-E hivpr subsample (133 ligands): ROC AUC 0.665.  This is an **EVALUATION** of the DockFlow+Vina pipeline on one target, NOT a claim of superior virtual-screening performance (single-target, exhaustiveness 1, no comparative baseline against other tools on the same panel).

Protocol: DUD-E hivpr actives + decoys (seeded subsample 40/536 actives, 93/35750 decoys of the full set), receptor 1HVR prepared by DockFlow (engine auto), grid box from the XK2 co-crystal ligand (strong assumption, 4.0 A padding), Vina python bindings, exhaustiveness 1, num_modes 1, seed 2026.  Ligand preparation: Meeko (embed 3D + MMFF minimise, seeded), input protonation states, no tautomer enumeration, Meeko-default salt stripping.  DockFlow example target; receptor + box reused from the redocking benchmark run (identical preparation).

- ligands scored: **133**
- ROC AUC: **0.665** (random baseline 0.499)
- EF@1%: 3.325 | EF@5%: 2.375 (random baseline 1.0)
- BEDROC (alpha 20): 0.6647 (random baseline 0.3001 - the BEDROC random baseline is NOT 0.5 at alpha 20; it depends on the actives fraction, so it is measured here over 10000 seeded random permutations of the same panel)

## Duplicate handling (item 10)

- DUD-E distributes already-deduplicated lists; this runner VERIFIES that (method: rdkit-canonical-smiles) instead of trusting it.  Duplicates are reported, never dropped - the published panel is exactly the seeded subsample.
- duplicate entries found: 4

## Interpretation

AUC > 0.5 means the Vina score ranks actives above decoys more often than chance; EF@x% is the actives-fold enrichment in the top x% of the ranking (1.0 = random).  All metrics are heuristic screening statistics, not binding free energies.  Caveats: single exhaustiveness-1 run, single target, no comparative baseline - the AUC is informative for 'does the pipeline enrich actives at all', not for 'is DockFlow better than tool X'.  Three targets is still a small panel; a publication-grade VS benchmark would use the full DUD-E (102 targets) - see docs/SCIENTIFIC_VALIDATION.md.

- ROC plot: `/home/z/my-project/dockflow-v030/benchmarks/enrichment/results/hivpr/roc.png`
- per-ligand scores: `/home/z/my-project/dockflow-v030/benchmarks/enrichment/results/hivpr/enrichment_results.csv`
- metrics JSON: `/home/z/my-project/dockflow-v030/benchmarks/enrichment/results/hivpr/enrichment_metrics.json`
