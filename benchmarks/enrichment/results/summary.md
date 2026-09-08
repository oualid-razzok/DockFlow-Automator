# DUD-E enrichment benchmark (target hivpr vs 1HVR)

Protocol: DUD-E hivpr actives + decoys (seeded subsample 40/536 actives, 93/35750 decoys of the full 536/35,750 set), docked with Vina (python bindings) against the 1HVR receptor, grid box from the XK2 co-crystal ligand (strong assumption), exhaustiveness 1, num_modes 1, seed 2026.

- ligands scored: **133**
- ROC AUC: **0.665**
- EF@1%: 3.325 | EF@5%: 2.375
- BEDROC (alpha 20): 0.6647

Interpretation: AUC > 0.5 means the Vina score ranks actives above decoys more often than chance; EF@x% is the actives-fold enrichment in the top x% of the ranking (1.0 = random). All metrics are heuristic screening statistics, not binding free energies.

- ROC plot: `/home/z/my-project/DockFlow-Automator/benchmarks/enrichment/results/roc.png`
- per-ligand scores: `/home/z/my-project/DockFlow-Automator/benchmarks/enrichment/results/enrichment_results.csv`
- metrics JSON: `/home/z/my-project/DockFlow-Automator/benchmarks/enrichment/results/enrichment_metrics.json`
