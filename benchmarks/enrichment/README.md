# Enrichment benchmark: actives vs decoys (audit item 22)

Quantifies the "decoys score worse" claim with a standardised
active/decoy panel instead of arbitrary aspirin/caffeine controls.

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
