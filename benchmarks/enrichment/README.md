# Enrichment benchmark: actives vs decoys (audit item 22)

Quantifies the "decoys score worse" claim with a standardised
active/decoy panel instead of arbitrary aspirin/caffeine controls.

## Status: RUN on the real DUD-E hivpr panel (2026-09)

`run_enrichment.py` docked a seeded subsample (40/536 actives,
93/35,750 decoys - single-CPU budget; ratio recorded in the summary)
against the 1HVR benchmark receptor with the XK2-derived box,
exhaustiveness 1, num_modes 1, seed 2026.  Results in `results/`:

* **ROC AUC 0.665** | EF@1 % 3.33 | EF@5 % 2.38 | BEDROC (α = 20) 0.665
  (133 ligands scored; per-ligand scores in `results/enrichment_results.csv`,
  ROC plot `results/roc.png`, machine-readable `results/enrichment_metrics.json`);
* interpretation: AUC > 0.5 means the Vina score ranks actives above
  decoys more often than chance.  These are **heuristic screening
  statistics, not binding free energies**;
* re-run / extend: `python benchmarks/enrichment/run_enrichment.py`
  (resumable; `--report-only` recomputes the metrics).

## Panel: DUD-E HIV-1 protease

The recommended panel is the **DUD-E** set for the example's target
(HIV-1 protease, `hivpr`): 39 actives + ~1 470 decoys with matched
physicochemistry.  Download (requires the DUD-E mirror; no registration):

```bash
python benchmarks/enrichment/download_dude.py --target hivpr --out benchmarks/enrichment/data
```

The script fetches `actives_final.ism` / `decoys_final.ism` (SMILES)
from `https://dude.docking.org/targets/<target>/`, prepares every
molecule with `dockflow prep ligand --library`, and writes a
`panel.sdf` plus an `actives.txt` id list ready for docking.

## Scoring workflow

```bash
# 1. dock the panel against the prepared 1HVR receptor (see example config)
dockflow dock --receptor runs/<run>/prepared/receptor.pdbqt \
    --ligands 'benchmarks/enrichment/data/prepared/*.pdbqt' \
    --config runs/<run>/gridbox.txt --out-dir enrichment_docking \
    --exhaustiveness 8 --seed 2026

# 2. enrichment metrics + ROC PNG
dockflow enrich --summary enrichment_docking/summary.csv \
    --actives benchmarks/enrichment/data/actives.txt \
    --roc-png enrichment_docking/roc.png --json
```

Metrics: **ROC AUC**, **EF@1 %**, **EF@5 %**, **BEDROC** (α = 20; the
random baseline is α- and actives-fraction-dependent, see
`dockflow_core/enrichment.py`).  `docking_score_efficiency` is a
Vina-score proxy, not experimental ligand efficiency.

## Aspirin/caffeine remain only as smoke controls

The 1HVR example config keeps aspirin/caffeine as **quick smoke-test
controls** (do they dock without errors and score worse than the
co-crystal ligand), not as enrichment evidence - with two decoys no
enrichment statistic is meaningful.
