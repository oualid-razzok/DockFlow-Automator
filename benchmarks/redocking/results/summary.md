# Redocking benchmark results

- complexes attempted: 24
- completed: 24 (0 failed - failures are recorded, not dropped)

## Pose recovery (formal definition: rank-1 pose, symmetry-aware heavy-atom RMSD <= 2.0 A)

- **rank-1 success: 21/24** = 87.5% (95% CI: 69.0%-95.7%)
- top-3 success (best RMSD among the 3 best-scored poses): **23/24** = 95.8% (95% CI: 79.8%-99.3%)
- best-any-rank success (v0.3.x legacy definition, kept for continuity): **23/24** = 95.8% (95% CI: 79.8%-99.3%)

Definition: 'pose recovery' = the best-SCORING Vina pose (rank 1) whose symmetry-aware heavy-atom RMSD to the co-crystallized ligand is <= 2.0 A.  The 2.0 A threshold is a community convention (Trott & Olson 2010; CSAR / D3R benchmark practice), NOT a thermodynamically validated cutoff - see docs/interaction_criteria.md.  Alternative definitions (best-RMSD-among-top-N) yield different success rates, which is exactly why all of them are published above.

## RMSD distribution (best pose per complex, A)

- min 0.02 | p25 0.26 | median 0.62 | p75 0.86 | max 3.03 | std 0.62 (n=24)
- histogram: `rmsd_distribution.png`
- mean 0.67 A | median 0.62 A

- mean runtime per complex: 58.5 s
- exhaustiveness: 8 (fixed seed 2026, success rates reported with Wilson score 95% confidence intervals; formula documented in run_benchmark.py)

## Success by definition and threshold (item 15)

| pose-recovery definition | <= 1.0 A | <= 2.0 A | <= 3.0 A |
|---|---|---|---|
| rank-1 (formal) | 15 | 21 | 22 |
| top-3 | 20 | 23 | 23 |
| best-any-rank (legacy) | 20 | 23 | 23 |

Top-N success curve (success rate when any of the first N scored poses is within 2.0 A): `topn_success.png`.  A flat curve for a complex means no good pose was sampled at all; a rising curve means a good pose exists but is outranked.

## Per-target grouping (item 2)

| target group | complexes | successes (best-any-rank, 2.0 A) | rate (95% CI) |
|---|---|---|---|
| HIV-1 protease | 6 | 5 | 83.3% (95% CI: 43.6%-97.0%) |
| other targets | 18 | 18 | 100.0% (95% CI: 82.4%-100.0%) |

## Failed complexes, legacy definition (best-any-rank, item 2)

Every failed complex had a strong ligand-derived grid box (see the per-complex table), so these failures are docking-side (sampling/scoring of flexible ligands), not grid-box or preparation-side.  Hypothesized causes below are the authors' reading of the provenance columns, not proven mechanisms.

| pdb | ligand | target | best RMSD (A) | rank-1 RMSD (A) | rotatable bonds | hypothesized cause |
|---|---|---|---|---|---|---|
| 1HTG | G37 | HIV-1 protease (penicillin-derived inhibitor) | 3.03 | 4.35 | 17 | penicillin-derived peptidomimetic, the most flexible ligand in the set (17 rotatable bonds): no pose within 3.0 A of the crystal pose among all 9 modes - a sampling limitation, not a ranking one (the top-N curve stays flat) |

## Failed complexes, formal definition (rank-1, item 2/15)

The formal pose-recovery definition is stricter: the best-SCORING pose must be within 2.0 A.  Distinguishing the two failure modes is the point of the top-N curve: a good pose that exists but ranks 3rd is a SCORING failure; no good pose at all is a SAMPLING failure.

| pdb | ligand | rank-1 RMSD (A) | best RMSD (A) at rank | failure mode |
|---|---|---|---|---|
| 1DWD | MID | 3.56 | 1.36 @ rank 3 | ranking (good pose sampled but outranked) |
| 1HPV | 478 | 2.67 | 0.69 @ rank 3 | ranking (good pose sampled but outranked) |
| 1HTG | G37 | 4.35 | 3.03 @ rank 4 | sampling (no pose within 2.0 A at any rank) |

## Per-complex results (with preparation + grid provenance, item 1)

Full provenance (grid-box center/size, docking parameters, RMSD method) is in `redocking_results.csv`; this table shows the columns a reader needs to judge whether a failure is prep- or docking-related.

| pdb | ligand | target | best affinity | rank-1 RMSD (A) | best RMSD (A) | success | gridbox source | assumption | engine | runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1DWD | MID | HIV-1 protease | -10.30 | 3.56 | 1.36 | yes | pocket:MID | strong | openbabel | 85.4 |
| 1F0R | 815 | human coagulation factor Xa (RPR208815) | -10.25 | 0.30 | 0.30 | yes | pocket:815 | strong | openbabel | 35.5 |
| 1SNC | THP | staphylococcal nuclease (thymidine 3'5'-diphosphate) | -6.92 | 1.60 | 0.99 | yes | pocket:THP | strong | openbabel | 33.7 |
| 1STP | BTN | streptavidin (biotin) | -7.43 | 0.62 | 0.62 | yes | pocket:BTN | strong | openbabel | 14.0 |
| 1CBX | BZS | carboxypeptidase A (benzylsuccinate) | -7.54 | 0.35 | 0.35 | yes | pocket:BZS | strong | openbabel | 12.0 |
| 1TNG | AMC | serine protease (benzamidine derivative) | -4.74 | 0.46 | 0.09 | yes | pocket:AMC | strong | openbabel | 5.5 |
| 1BTY | BEN | beta-trypsin (benzamidine) | -5.55 | 0.14 | 0.11 | yes | pocket:BEN | strong | openbabel | 5.6 |
| 1MQ6 | XLD | 3-chloro-pyridyl inhibitor complex | -9.52 | 1.62 | 0.77 | yes | pocket:XLD | strong | openbabel | 58.3 |
| 1JAP | HOA | pro-leu-gly-hydroxylamine complex | -2.38 | 0.02 | 0.02 | yes | pocket:HOA | strong | openbabel | 4.5 |
| 1YGC | 905 | factor VIIa (small molecule inhibitor) | -9.78 | 0.65 | 0.65 | yes | pocket:905 | strong | openbabel | 82.6 |
| 1SL3 | 170 | thrombin (P1 heteroaryl inhibitor) | -10.18 | 0.97 | 0.97 | yes | pocket:170 | strong | openbabel | 79.3 |
| 1HSG | MK1 | HIV-1 protease (L-735,524 saquinavir analog) | -11.34 | 0.42 | 0.42 | yes | pocket:MK1 | strong | openbabel | 128.3 |
| 1HPV | 478 | HIV-1 protease (VX-478 amprenavir) | -8.83 | 2.67 | 0.69 | yes | pocket:478 | strong | openbabel | 64.9 |
| 1MUI | AB1 | HIV-1 protease (lopinavir) | -10.54 | 0.47 | 0.47 | yes | pocket:AB1 | strong | openbabel | 165.2 |
| 1AZM | AZM | human carbonic anhydrase II (acetazolamide) | -6.08 | 1.13 | 0.62 | yes | pocket:AZM | strong | openbabel | 6.8 |
| 1E66 | HUX | acetylcholinesterase (huprine X) | -11.40 | 0.26 | 0.26 | yes | pocket:HUX | strong | openbabel | 8.8 |
| 3ERT | OHT | human estrogen receptor alpha (4-hydroxytamoxifen) | -9.71 | 1.00 | 0.83 | yes | pocket:OHT | strong | openbabel | 37.5 |
| 1N8Q | DHB | lipoxygenase (protocatechuic acid) | -5.56 | 0.32 | 0.22 | yes | pocket:DHB | strong | openbabel | 7.4 |
| 1NVQ | UCN | checkpoint kinase Chk1 (UCN-01) | -12.51 | 1.13 | 1.13 | yes | pocket:UCN | strong | openbabel | 16.9 |
| 1Q41 | IXM | GSK-3 beta (indirubin-3'-monoxime) | -10.56 | 0.25 | 0.25 | yes | pocket:IXM(A451) | strong | openbabel | 7.1 |
| 1XLZ | FIL | human phosphodiesterase 4B (filaminast) | -8.19 | 1.29 | 1.01 | yes | pocket:FIL(A1003) | strong | openbabel | 16.0 |
| 1N46 | PFA | thyroid hormone receptor beta | -11.53 | 0.18 | 0.18 | yes | pocket:PFA(A462) | strong | openbabel | 17.9 |
| 1HVR | XK2 | HIV-1 protease (cyclic urea inhibitor) | -11.18 | 0.97 | 0.81 | yes | pocket:XK2(A263) | strong | openbabel | 242.6 |
| 1HTG | G37 | HIV-1 protease (penicillin-derived inhibitor) | -10.56 | 4.35 | 3.03 | NO | pocket:G37(A300) | strong | openbabel | 269.0 |
