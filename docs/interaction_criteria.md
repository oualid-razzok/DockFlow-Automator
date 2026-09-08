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

## Redocking RMSD conventions (audit item 18, peer item 17)

`crystal_rmsd` in `summary.csv`, `interactions.json` and `report.md` is:

* **heavy-atom** (hydrogens excluded on both sides),
* **symmetry-aware** (peer item 17): with RDKit available the RMSD is
  `GetBestRMS` — the exact minimum over all automorphisms of the
  molecular graph (carboxylate flips, ring-atom permutations, terminal
  group swaps), which is the standard symmetry-correct ligand RMSD;
  without RDKit a greedy element-matching Kabsch fallback is used.  The
  method actually used is recorded per pose
  (`interactions.json` → `poses[*].crystal_rmsd_method`), per run
  (`manifest.analysis.crystal_rmsd_method`) and per benchmark row
  (`redocking_results.csv` → `crystal_rmsd_method`),
* **atom-order independent** (final-pass item 5): the pose molecule is
  built with its own proximity-bonded graph, so the RMSD does not depend
  on the order in which the pose file or the reference PDB lists atoms.
  This was not true before v1.0.0-rc1 — a geometrically seeded
  correspondence could lock in a topologically wrong assignment and
  over-penalise symmetric ligands by up to ~0.5 Å.  The fix was found by
  the independent spyrmsd cross-check
  (`benchmarks/redocking/rmsd_crosscheck.py`; 193 poses, max |Δ| = 0.0 Å
  after the fix) and is guarded by analytic test cases
  (`tests/test_analyzer.py::test_symmetry_rmsd_known_cases`),
* computed after optimal **superposition** (rotation + translation
  removed; reflections are not considered),
* reported against the co-crystallized ligand pose extracted from the
  target PDB (highest-occupancy altLoc resolved),
* undefined (`null`) when the docked ligand and the reference ligand have
  different heavy-atom element sets (e.g. decoys) — never an error.

## Pose recovery: the formal definition (final pass, item 3)

**"Pose recovery"** = the best-SCORING Vina pose (rank 1, the pose Vina
assigns the lowest predicted affinity) whose symmetry-aware heavy-atom
RMSD to the co-crystallized ligand is **≤ 2.0 Å**.

Two things a reader must know about this definition:

1. **The 2.0 Å threshold is a community convention, not a validated
   cutoff.**  It follows common redocking-benchmark practice since Trott
   & Olson (2010) and is the working convention of the CSAR and D3R
   benchmarking exercises.  It is NOT a thermodynamically validated
   boundary — a 1.9 Å pose is not "twice as good" as a 3.8 Å pose, and
   for some targets the crystal pose itself is one of several
   well-populated conformers.
2. **"Best-scoring" means Vina rank 1.**  Alternative definitions yield
   different success rates, which is why the benchmark publishes ALL of
   them (rank-1, top-3 best-RMSD, and the legacy best-any-rank, each at
   1.0/2.0/3.0 Å thresholds, plus the top-N success curve):

   * `rank1_rmsd_a` — RMSD of the rank-1 pose (the formal definition);
   * `top3_rmsd_a` — best RMSD among the 3 best-scored poses ("a good
     pose exists and Vina ranks it within the top 3");
   * `best_crystal_rmsd_a` — best RMSD among all returned poses (the
     v0.3.x legacy definition; sampling capability, not ranking).

   On the 24-complex benchmark (v1.0.0-rc1 numbers): rank-1 21/24
   (87.5%), top-3 23/24 (95.8%), best-any-rank 23/24 (95.8%) — see
   `benchmarks/redocking/results/summary.md` for the confidence
   intervals.  The distinction matters: a complex that fails rank-1 but
   passes top-3 is a *scoring* failure (the pose was found but
   outranked); a complex that fails at every rank is a *sampling*
   failure.

## Threshold classification: conventions vs heuristics (peer item 5; none experimentally validated)

Every threshold DockFlow applies falls into one of two classes.  The
report labels heuristic outputs explicitly so no number is mistaken for a validated result (the labels are heuristic, none experimentally confirmed).

| Threshold | Value | Class | Where used |
|---|---|---|---|
| Redocking success | ≤ 2.0 Å | **Convention** (community redocking practice; not validated for any specific target) | report banner, benchmark summary |
| Pose clustering | 2.0 Å Kabsch RMSD | **Heuristic** (no thermodynamic meaning; cluster membership is geometry only) | report "Heuristic pose clustering" |
| Docking score efficiency | affinity / heavy atoms | **Heuristic proxy** (Vina-score-derived, not experimental ligand efficiency) | report results table |
| Geometric H-bond contact | ≤ 3.5 Å | **Geometric convention** (distance only, no angle) | contacts CSV/JSON |
| Hydrophobic contact | ≤ 4.5 Å | **Geometric convention** | contacts CSV/JSON |
| Ionic contact | ≤ 4.5 Å | **Geometric convention** (residue-name qualified) | contacts CSV/JSON |
| Metal contact | ≤ 3.0 Å | **Geometric convention** | contacts CSV/JSON |
| Disulfide detection | S–S < 2.5 Å | **Geometric convention** (standard S–S bond length range) | manifest, warnings |
| Metal coordination | ≤ 2.8 Å | **Heuristic** (coordination number only; geometry label is a count heuristic) | manifest, warnings |
| Reduced Cys flagging | free SG within 8 Å of free SG | **Heuristic** (no redox context modelled) | manifest, warnings |
| Grid box ligand padding | 4 Å default | **Convention** (covers typical pose spread; no target validation) | grid box derivation |
| Box volume sanity | > 8 000 Å³ warn / > 27 000 Å³ error | **Heuristic guardrail** (typical binding-site scale) | grid box validation |

None of these values are experimentally validated for a specific target (they are conventions and heuristics, none
experimentally confirmed).  The only numbers in a DockFlow run that come from
outside this codebase are the Vina/GNINA scoring-function parameters,
and those are *empirical scoring functions*, not experimental
measurements.  When a threshold matters for a decision (publication,
go/no-go on a series), validate it against reference ligands for that
target — e.g. `benchmarks/redocking/` shows the pipeline's own
pose-recovery distribution on 24 public complexes.
