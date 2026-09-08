# DockFlow-Automator run report

- run id: **1HVR**
- target: **1HVR** (pdb), RATIONAL DESIGN OF POTENT, BIOAVAILABLE, NONPEPTIDE CYCLIC UREAS AS HIV PROTEASE INHIBITORS
- receptor: `examples/results/1HVR/prepared/receptor.pdbqt` (1846 atoms, engine: openbabel)
- grid box: center (-8.71, 15.54, 27.9), size (17.97, 15.46, 25.1) A (source: pocket:XK2 chain A residue 263, padding 4.0 A, assumption: strong)
- docking backend: python

> Comparability: results from this run are NOT directly comparable to receptors prepared with rdkit, none (hydrogen/charge differences); see manifest.receptor.decisions.comparable_to.

> Protonation: Catalytic residues detected near the binding site; protonation states may be wrong without pKa analysis. Hydrogens were added by toolkit template (pH 7.4 nominal), NOT by a pKa calculation (pH-sensitive residues near the pocket: GLU A21, ASP A25, ASP A29, ASP A30, GLU A34, GLU A35, GLU B21, ASP B25, ASP B29, ASP B30, GLU B34, GLU B35). Run PROPKA / H++ and pass the protonated PDB as the receptor input for rigorous work.

> **Redocking pose recovery: best pose RMSD = 0.84 A for xk2_redock** - PASS

RMSD is symmetry-aware heavy-atom RMSD to the co-crystallized ligand pose extracted from the target structure Method: rdkit-getbestrms-graph-automorphism.
The 2.0 A pass mark is a **heuristic success threshold** (community convention, not a validated cutoff for this target); see docs/interaction_criteria.md.

## Results

| ligand | best affinity (kcal/mol) | best crystal RMSD (A) | poses | docking score efficiency (heuristic proxy) | runtime (s) |
|---|---|---|---|---|---|
| xk2_redock | -11.13 | 0.84 | 9 | -0.242 | 463.8 |
| aspirin_decoy | -5.47 | n/a | 9 | -0.421 | 7.8 |
| caffeine_decoy | -5.43 | n/a | 9 | -0.388 | 3.6 |

_docking score efficiency (heuristic proxy) = affinity / heavy-atom count: a Vina-score-derived proxy, **not** experimental ligand efficiency._

## Heuristic pose clustering (xk2_redock)

Heuristic clustering: Kabsch RMSD, 2.0 A threshold; no thermodynamic meaning (cluster membership says nothing about binding free energy).

9 poses in **4 cluster(s)** (multiple distinct binding modes)

| cluster | poses | mean dG (kcal/mol) | spread (A) | representative pose |
|---|---|---|---|---|
| 1 | 1, 2, 3, 5, 6 | -10.22 | 2.25 | 1 |
| 2 | 4 | -9.81 | 0.00 | 4 |
| 3 | 7 | -9.12 | 0.00 | 7 |
| 4 | 8, 9 | -8.91 | 1.92 | 8 |

## Residue contact hotspots (best ligand, top pose)

Contact types are **geometric criteria** (distance + atom type only, no angle criterion); see docs/interaction_criteria.md.

| chain | residue | geom. H-bond contacts (<=3.5 A) | hydrophobic | ionic | metal | closest A |
|---|---|---|---|---|---|---|
| A | ILE50 | 1 | 16 | 0 | 0 | 2.37 |
| A | ILE47 | 0 | 17 | 0 | 0 | 3.56 |
| B | ILE47 | 0 | 11 | 0 | 0 | 3.81 |
| B | GLY49 | 0 | 10 | 0 | 0 | 3.47 |
| A | ALA28 | 0 | 9 | 0 | 0 | 3.34 |
| B | ILE50 | 0 | 9 | 0 | 0 | 3.69 |
| A | GLY49 | 0 | 7 | 0 | 0 | 3.78 |
| A | ILE84 | 0 | 7 | 0 | 0 | 3.89 |
| B | ALA28 | 1 | 5 | 0 | 0 | 3.44 |
| A | ASP30 | 0 | 6 | 0 | 0 | 3.71 |
| B | ASP30 | 0 | 6 | 0 | 0 | 3.83 |
| B | ASP25 | 4 | 0 | 0 | 0 | 2.12 |
| B | VAL32 | 0 | 4 | 0 | 0 | 3.57 |
| A | VAL82 | 0 | 3 | 0 | 0 | 3.75 |
| B | ILE84 | 0 | 3 | 0 | 0 | 3.75 |

## Warnings

- receptor: removed hetero species: 18 CSO, 46 XK2. These may be structurally or catalytically required (metals/cofactors); keep them with receptor.keep_resnames or keep_hetero if they belong in the binding site
- receptor: Receptor has 2 missing residues (chain A: 67; chain B: 67). Vina cannot model loops; consider Modeller/SWISS-MODEL before docking. (Heuristic numbering-gap detection; see manifest.receptor.decisions.missing_residues.)

## Assumptions this run made

- receptor treated as rigid: yes
- protonation: as deposited (PDB); polar hydrogens added by toolkit template (openbabel engine, pH 7.4 nominal; NOT a pKa calculation) - catalytic residues may be mis-protonated; for rigorous work, run PROPKA and pass the protonated PDB as input
- charge model: gasteiger (Gasteiger partial charges)
- waters: removed (0 water oxygen atoms)
- metals/cofactors: removed {'CSO': 18, 'XK2': 46} / kept none
- grid box: derived from pocket:XK2(A263), padding 4.0 A (assumption: strong)
- scoring: vina (AutoDock Vina 1.2.7, no rescoring)
- ligand tautomers/protonation: see manifest ligands[*].prep (input state by default; no enumeration unless configured)

See README section *Scientific limitations and assumptions* for the general limitations behind each of these decisions, and docs/preparation_assumptions.md for the full preparation assumption table (what DockFlow does / does NOT do / what the researcher should do for publication-grade work).

## Files

- manifest: `examples/results/1HVR/manifest.json`
- environment fingerprint: `examples/results/1HVR/environment.json`
- docking summary: `examples/results/1HVR/docking/summary.csv`
- interactions: `examples/results/1HVR/analysis/interactions.json`
- visualization: `examples/results/1HVR/visualization`

_Generated by DockFlow-Automator v0.3.1_
