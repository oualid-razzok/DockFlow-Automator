# DockFlow vs MGLTools preparation comparison

Both arms: AutoDock Vina (python bindings), the SAME search space
(redocking benchmark gridbox), seed 2026, exhaustiveness
8, num_modes 9.  The only difference is the
preparation engine - every delta below is attributable to
preparation chemistry (hydrogen model, charge assignment), not to
the docking engine.

- complexes attempted: 24
- MGLTools arm completed: 21 (failures recorded with reasons, not dropped)
- MGLTools arm success (RMSD <= 2.0 A): **18/21** = 85.7% (95% CI: 65.4%-95.0%)
- DockFlow arm success: **20/21** = 95.2% (95% CI: 77.3%-99.2%)
- **denominator (explicit, item 4): the success rates above are computed on the 21-complex MGLTools-completed subset - 3 MGLTools crashes are excluded.  This is NOT a head-to-head 24-vs-24 comparison; the full 24-complex DockFlow numbers are in benchmarks/redocking/results/.
- mean RMSD delta (MGLTools - DockFlow): 0.253 A
- mean affinity delta (MGLTools - DockFlow): -0.087 kcal/mol
- mean receptor atom-count delta: 0.0

## Per-complex deltas

| pdb | ligand | RMSD DockFlow (A) | RMSD MGLTools (A) | delta | affinity DockFlow | affinity MGLTools | delta | rec atoms DF/MGL | lig atoms DF/MGL |
|---|---|---|---|---|---|---|---|---|---|
| 1HVR | XK2 | 0.81 | 3.19 | 2.37 | -11.18 | -12.48 | -1.30 | 1500/1500 | 44/46 |
| 1HSG | MK1 | 0.42 | 0.61 | 0.19 | -11.34 | -11.38 | -0.05 | 1514/1514 | 45/45 |
| 1HTG | G37 | 3.03 | 2.10 | -0.93 | -10.56 | -10.26 | 0.31 | 1516/1516 | 55/55 |
| 1HPV | 478 | 0.69 | 1.33 | 0.64 | -8.83 | -8.88 | -0.04 | 1516/1516 | 35/35 |
| 1MUI | AB1 | 0.47 | 0.84 | 0.37 | -10.54 | -11.03 | -0.48 | 1516/1516 | 46/46 |
| 1DWD | MID | 1.36 | 1.37 | 0.01 | -10.30 | -9.83 | 0.47 | 2358/2358 | 37/37 |
| 1AZM | AZM | 0.62 | 0.57 | -0.06 | -6.08 | -6.27 | -0.19 | 2019/2019 | 13/13 |
| 1E66 | HUX | 0.26 | FAILED (receptor_prep) | - | -11.40 | - | - | - | - | root cause: IndexError: list index out of range; [full traceback](full traceback: failures/1E66_mgltools_receptor_prep_traceback.txt) |
| 1F0R | 815 | 0.30 | 0.33 | 0.03 | -10.25 | -10.34 | -0.09 | 2227/2227 | 31/31 |
| 1Q41 | IXM | 0.25 | 0.38 | 0.13 | -10.56 | -10.22 | 0.34 | 5428/5428 | 21/21 |
| 1SNC | THP | 0.99 | 1.09 | 0.09 | -6.92 | -7.40 | -0.48 | 1082/1082 | 25/25 |
| 3ERT | OHT | 0.83 | 1.05 | 0.22 | -9.71 | -9.54 | 0.16 | 1932/1932 | 29/29 |
| 1XLZ | FIL | 1.01 | 1.00 | -0.02 | -8.19 | -8.26 | -0.07 | 5202/5202 | 21/21 |
| 1STP | BTN | 0.62 | 0.60 | -0.02 | -7.43 | -7.24 | 0.19 | 901/901 | 16/16 |
| 1CBX | BZS | 0.35 | 0.17 | -0.18 | -7.54 | -8.02 | -0.48 | 2437/2437 | 15/15 |
| 1TNG | AMC | 0.09 | 0.16 | 0.07 | -4.74 | -4.74 | 0.00 | 1629/1629 | 8/8 |
| 1BTY | BEN | 0.11 | 0.05 | -0.05 | -5.55 | -5.43 | 0.12 | 1629/1629 | 9/9 |
| 1MQ6 | XLD | 0.77 | FAILED (ligand_prep) | - | -9.52 | - | - | - | - | root cause: IndexError: string index out of range; [full traceback](full traceback: failures/1MQ6_mgltools_ligand_prep_traceback.txt) |
| 1N8Q | DHB | 0.22 | 0.25 | 0.03 | -5.56 | -5.66 | -0.10 | 6778/6778 | 11/11 |
| 1NVQ | UCN | 1.13 | 2.56 | 1.44 | -12.51 | -11.89 | 0.61 | 2169/2169 | 36/36 |
| 1JAP | HOA | 0.02 | 0.03 | 0.00 | -2.38 | -2.38 | 0.00 | 1260/1260 | 2/2 |
| 1YGC | 905 | 0.65 | 1.40 | 0.74 | -9.78 | -10.07 | -0.29 | 2368/2368 | 38/38 |
| 1N46 | PFA | 0.18 | 0.39 | 0.22 | -11.53 | -12.01 | -0.48 | 3872/3872 | 27/27 |
| 1SL3 | 170 | 0.97 | FAILED (ligand_prep) | - | -10.18 | - | - | - | - | root cause: IndexError: string index out of range; [full traceback](full traceback: failures/1SL3_mgltools_ligand_prep_traceback.txt) |

## Interpretation

MGLTools and DockFlow (openbabel engine) differ in hydrogen
placement and charge assignment, so the same Vina search lands on
slightly different poses and scores.  A positive RMSD delta means
the MGLTools arm recovered the crystal pose better; negative means
DockFlow did.  Neither arm is a ground truth: both are compared
against the crystal pose.

The MGLTools-arm failures are **genuine MGLTools 1.5.7 crashes** on these inputs (full tracebacks in ``results/failures/``, linked from the table and the CSV); the DockFlow arm processed the same files successfully.  They are recorded, not dropped, so the two arms' success rates share the MGLTools-completed denominator stated above - not a head-to-head 24-vs-24.

Systematic differences to keep in mind when reading the table:

* **Atom counts**: the receptor atom-count delta counts non-hydrogen
  PDBQT atoms (virtual glue atoms excluded); MGLTools' ``-U nphs``
  removes non-polar hydrogens (united-atom model, polar H kept),
  DockFlow's openbabel engine also keeps polar hydrogens only - small
  count differences reflect H-bonding assignment differences, not
  lost atoms.  The ligand column counts real heavy atoms; polar-H
  differences are in ``mgl_ligand_polar_h``.
* **Ligand torsions**: ``prepare_ligand4`` and Meeko can disagree
  on which bonds are rotatable, which changes the search-space size.
* **Macrocycles**: Meeko opens macrocyclic rings for flexible-
  macrocycle docking by replacing the two ring-closure atoms with
  G0/CG0 glue pseudo-atoms - for such ligands (1HVR/XK2 in this set)
  the DockFlow-arm explicit heavy count is 2 lower than MGLTools'
  (which docks the closed ring); no atoms are lost.
* **Charges**: both arms use Gasteiger charges, but hydrogen
  placement changes how those charges are distributed.

