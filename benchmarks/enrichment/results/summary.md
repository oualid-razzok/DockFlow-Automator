# DUD-E enrichment evaluation - multi-target summary (item 12)

| target | receptor (PDB) | n actives | n decoys | ROC AUC | EF@1% | EF@5% | BEDROC (a=20) |
|---|---|---|---|---|---|---|---|
| hivpr | 1HVR | 40 | 93 | 0.665 | 3.325 | 2.375 | 0.6647 |
| aa2ar | 3EML | 40 | 93 | 0.593 | 1.663 | 1.9 | 0.5486 |
| parp1 | 2RD6 | 40 | 93 | 0.729 | 3.325 | 1.425 | 0.5157 |

All targets share the identical protocol (see each target's summary.md): seeded 40/93 subsample, Meeko ligand prep, ligand-derived grid box (4.0 A padding), Vina exhaustiveness 1, num_modes 1, seed 2026.  ROC AUC random baseline 0.500; EF random baseline 1.0; the BEDROC random baseline depends on the actives fraction (recorded per target in enrichment_metrics.json).

**Limitation (explicit):** three targets is still a small panel and every number above is a single exhaustiveness-1 seeded run.  A publication-grade virtual-screening benchmark would evaluate the full DUD-E (102 targets) with multiple seeds - that is future work (see ROADMAP_POST_V1.md), stated here so the table is not over-read.

