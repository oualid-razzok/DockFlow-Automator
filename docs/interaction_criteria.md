# Interaction (contact) criteria

DockFlow's ligand–receptor "interaction" analysis is **geometric /
contact-based**: it classifies atom pairs by *distance and AutoDock atom
type only*. It is not a rigorous hydrogen-bond, ionic-bond or
metal-coordination analysis — no angle or orbital-geometry criteria are
applied. Every user-facing surface therefore reports **"geometric H-bond
contact"** (etc.), never "H-bond".

This page documents the exact criteria so results are interpretable and
reproducible, and states the limitations of the model.

## Criteria table

| Classification | Ligand atom (PDBQT type) | Receptor atom (PDBQT type) | Distance cutoff | Additional requirement |
|---|---|---|---|---|
| Geometric H-bond contact | `HD` (polar hydrogen) | `NA OA SA OS N O S` | ≤ 3.5 Å | none |
| Geometric H-bond contact | `NA OA SA OS N O S` | `HD` | ≤ 3.5 Å | none |
| Hydrophobic contact | `C A` | `C A` | ≤ 4.5 Å | none |
| Ionic contact | ligand cation (`N NA NS`) | ASP/GLU `OD*`/`OE*` carboxylate O | ≤ 4.5 Å | receptor residue name |
| Ionic contact | ligand anion (`OA OS O SA S`) | ARG/LYS `NH*`/`NZ` protonated N | ≤ 4.5 Å | receptor residue name |
| Metal contact | ligand N/O/S types | receptor metal (`Mg Ca Mn Fe Zn NI CU CO`) | ≤ 3.0 Å | none |

Implementation: `dockflow_core/analyzer.py` (`classify_pair`,
`HBOND_CUTOFF`, `HYDROPHOBIC_CUTOFF`, `IONIC_CUTOFF`, `METAL_CUTOFF`).

Constants of the model:

* Contact detection cutoff (overall): 5.0 Å (`CONTACT_CUTOFF`,
  configurable via `analysis.contacts_cutoff`).
* Distances are computed between **PDBQT atom coordinates**; the receptor
  keeps only polar hydrogens (`HD`) because non-polar hydrogens are merged
  into their heavy atoms during preparation (as in
  `prepare_receptor4.py -U nphs`).

## What this model can and cannot tell you

**Can**: which residues line the binding site; which pairs are close
enough for the interaction type to be *possible*; a residue-level ranking
of contact density; reproducible, deterministic contact fingerprints for
comparing poses and ligands.

**Cannot**:

* Distinguish a genuine, geometrically well-formed hydrogen bond from a
  merely short donor–acceptor contact (no D–H···A angle criterion is
  applied; a 3.4 Å contact at a 90° D–H···A angle is counted the same as
  one at 170°).
* Detect directionality of ionic interactions or metal coordination
  geometry.
* Say anything about interaction *strength* — contact counts are not
  energies.
* Detect water-mediated contacts (waters are removed during receptor
  preparation by default).

If angle-qualified interaction perception is required, post-process the
docked PDBQT with a dedicated tool (e.g. PLIP, Arpeggio) — DockFlow's CSV
and JSON contact tables are designed to be importable for such workflows.

## Where these classifications appear

* `analysis/*_poseN_contacts.csv` — one row per contact pair (columns:
  `kind`, `distance`, atom names/types; `kind` uses the short tags
  `hbond`/`hydrophobic`/`ionic`/`metal`, all meaning the geometric
  criteria above).
* `analysis/interactions.json` — per-pose contact lists and per-residue
  aggregates; keys are explicitly qualified
  (`num_geom_hbond_contacts`, `geom_hbond_contacts`, ...). The file
  validates against `docs/schemas/interactions.schema.json`.
* `report.md` — "Residue contact hotspots" table with qualified column
  headings and cutoffs in the header.
* CLI (`dockflow analyze`) — prints "geometric H-bond contacts" counts.

## Redocking RMSD conventions (audit item 18)

`crystal_rmsd` in `summary.csv`, `interactions.json` and `report.md` is:

* **heavy-atom** (hydrogens excluded on both sides),
* **symmetry-tolerant**: chemically equivalent atoms (same element) may
  be permuted — a 180° flip of a symmetric carboxylate or ring nitrogen
  permutation costs 0 Å, matching standard redocking-benchmark practice,
* computed after optimal **Kabsch superposition** (rotation + translation
  removed),
* reported against the co-crystallized ligand pose extracted from the
  target PDB (highest-occupancy altLoc resolved),
* undefined (`null`) when the docked ligand and the reference ligand have
  different heavy-atom element sets (e.g. decoys) — never an error.

Success threshold for pose recovery: best pose RMSD ≤ 2.0 Å (the
conventional redocking criterion; the report banner states it
explicitly).
