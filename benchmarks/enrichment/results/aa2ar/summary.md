# DUD-E enrichment evaluation (target aa2ar, subsample 133 ligands)

Enrichment evaluation on the DUD-E aa2ar subsample (133 ligands): ROC AUC 0.593.  This is an **EVALUATION** of the DockFlow+Vina pipeline on one target, NOT a claim of superior virtual-screening performance (single-target, exhaustiveness 1, no comparative baseline against other tools on the same panel).

Protocol: DUD-E aa2ar actives + decoys (seeded subsample 40/482 actives, 93/31550 decoys of the full set), receptor 3EML prepared by DockFlow (engine auto), grid box from the ZMA co-crystal ligand (strong assumption, 4.0 A padding), Vina python bindings, exhaustiveness 1, num_modes 1, seed 2026.  Ligand preparation: Meeko (embed 3D + MMFF minimise, seeded), input protonation states, no tautomer enumeration, Meeko-default salt stripping.  adenosine A2A receptor; reference structure 3EML (ZM241385 co-crystal).

- ligands scored: **133**
- ROC AUC: **0.593** (random baseline 0.500)
- EF@1%: 1.663 | EF@5%: 1.9 (random baseline 1.0)
- BEDROC (alpha 20): 0.5486 (random baseline 0.3008 - the BEDROC random baseline is NOT 0.5 at alpha 20; it depends on the actives fraction, so it is measured here over 10000 seeded random permutations of the same panel)

## Duplicate handling (item 10)

- DUD-E distributes already-deduplicated lists; this runner VERIFIES that (method: rdkit-canonical-smiles) instead of trusting it.  Duplicates are reported, never dropped - the published panel is exactly the seeded subsample.
- duplicate entries found: 0

## Interpretation

AUC > 0.5 means the Vina score ranks actives above decoys more often than chance; EF@x% is the actives-fold enrichment in the top x% of the ranking (1.0 = random).  All metrics are heuristic screening statistics, not binding free energies.  Caveats: single exhaustiveness-1 run, single target, no comparative baseline - the AUC is informative for 'does the pipeline enrich actives at all', not for 'is DockFlow better than tool X'.  Three targets is still a small panel; a publication-grade VS benchmark would use the full DUD-E (102 targets) - see docs/SCIENTIFIC_VALIDATION.md.

- ROC plot: `/home/z/my-project/dockflow-v030/benchmarks/enrichment/results/aa2ar/roc.png`
- per-ligand scores: `/home/z/my-project/dockflow-v030/benchmarks/enrichment/results/aa2ar/enrichment_results.csv`
- metrics JSON: `/home/z/my-project/dockflow-v030/benchmarks/enrichment/results/aa2ar/enrichment_metrics.json`
