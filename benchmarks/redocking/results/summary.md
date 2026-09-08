# Redocking benchmark results

- complexes attempted: 24
- completed: 24 (0 failed - failures are recorded, not dropped)
- successful (best pose RMSD <= 2.0 A, heuristic community convention): **18/24** (75%)
- mean best-pose RMSD: 1.32 A | median: 0.96 A
- mean runtime per complex: 58.5 s
- exhaustiveness: 8 (fixed seed 2026)

## Per-complex results

| pdb | ligand | target | best affinity | best RMSD (A) | success | runtime (s) | engine |
|---|---|---|---|---|---|---|---|
| 1DWD | MID | HIV-1 protease | -10.30 | 2.52 | NO | 85.4 | openbabel |
| 1F0R | 815 | human coagulation factor Xa (RPR208815) | -10.25 | 2.69 | NO | 35.5 | openbabel |
| 1SNC | THP | staphylococcal nuclease (thymidine 3'5'-diphosphate) | -6.92 | 0.99 | yes | 33.7 | openbabel |
| 1STP | BTN | streptavidin (biotin) | -7.43 | 0.62 | yes | 14.0 | openbabel |
| 1CBX | BZS | carboxypeptidase A (benzylsuccinate) | -7.54 | 0.35 | yes | 12.0 | openbabel |
| 1TNG | AMC | serine protease (benzamidine derivative) | -4.74 | 0.09 | yes | 5.5 | openbabel |
| 1BTY | BEN | beta-trypsin (benzamidine) | -5.55 | 0.11 | yes | 5.6 | openbabel |
| 1MQ6 | XLD | 3-chloro-pyridyl inhibitor complex | -9.52 | 2.86 | NO | 58.3 | openbabel |
| 1JAP | HOA | pro-leu-gly-hydroxylamine complex | -2.38 | 0.02 | yes | 4.5 | openbabel |
| 1YGC | 905 | factor VIIa (small molecule inhibitor) | -9.78 | 3.67 | NO | 82.6 | openbabel |
| 1SL3 | 170 | thrombin (P1 heteroaryl inhibitor) | -10.18 | 0.97 | yes | 79.3 | openbabel |
| 1HSG | MK1 | HIV-1 protease (L-735,524 saquinavir analog) | -11.34 | 4.01 | NO | 128.3 | openbabel |
| 1HPV | 478 | HIV-1 protease (VX-478 amprenavir) | -8.83 | 0.95 | yes | 64.9 | openbabel |
| 1MUI | AB1 | HIV-1 protease (lopinavir) | -10.54 | 1.47 | yes | 165.2 | openbabel |
| 1AZM | AZM | human carbonic anhydrase II (acetazolamide) | -6.08 | 0.62 | yes | 6.8 | openbabel |
| 1E66 | HUX | acetylcholinesterase (huprine X) | -11.40 | 0.26 | yes | 8.8 | openbabel |
| 3ERT | OHT | human estrogen receptor alpha (4-hydroxytamoxifen) | -9.71 | 1.88 | yes | 37.5 | openbabel |
| 1N8Q | DHB | lipoxygenase (protocatechuic acid) | -5.56 | 0.68 | yes | 7.4 | openbabel |
| 1NVQ | UCN | checkpoint kinase Chk1 (UCN-01) | -12.51 | 0.84 | yes | 16.9 | openbabel |
| 1Q41 | IXM | GSK-3 beta (indirubin-3'-monoxime) | -10.56 | 0.25 | yes | 7.1 | openbabel |
| 1XLZ | FIL | human phosphodiesterase 4B (filaminast) | -8.19 | 1.01 | yes | 16.0 | openbabel |
| 1N46 | PFA | thyroid hormone receptor beta | -11.53 | 0.18 | yes | 17.9 | openbabel |
| 1HVR | XK2 | HIV-1 protease (cyclic urea inhibitor) | -11.18 | 1.06 | yes | 242.6 | openbabel |
| 1HTG | G37 | HIV-1 protease (penicillin-derived inhibitor) | -10.56 | 3.64 | NO | 269.0 | openbabel |
