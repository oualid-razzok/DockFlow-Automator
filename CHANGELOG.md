# Changelog

## [0.3.0] - 2026-09-08

Validation-track release: the second-round work order (31 items).  The
headline change is scientific honesty - the validation suite that 0.2.1
*claimed* now actually *ran*, and every heuristic the pipeline applies
is labelled as one.

### Validation track (work order items 1-3) - RUN, not promised
- **24-complex redocking benchmark: complete.**  All 24 complexes
  docked (exhaustiveness 8, seed 2026): 24/24 completed, **18/24
  successes** (75%, best pose RMSD ≤ 2.0 Å), mean best-pose RMSD
  1.32 Å / median 0.96 Å; per-complex table, scatter plot and
  per-complex artifacts in `benchmarks/redocking/results/`.  Six fails
  are all flexible peptidomimetic HIV-PR inhibitors (recorded, kept,
  inspectable).  Runner is resumable/crash-safe with portable
  relative run dirs; two independent executions reproduced identical
  RMSDs (seed determinism).
- **MGLTools comparison: run for real.**  MGLTools 1.5.7 (bioconda
  package with bundled Python 2.7, relocated in-place via
  `scripts/fix_mgltools_prefix.py`) prepared the conventional arm on
  the same inputs, box and seed: 21/24 completed, **3 genuine MGLTools
  1.5.7 crashes** (recorded with tracebacks, never dropped); pose
  recovery DockFlow 16/21 vs MGLTools 14/21, mean RMSD delta
  +0.45 Å, receptor atom counts identical on every complex.
- **DUD-E enrichment: run on the real hivpr panel.**  133 ligands
  (40/536 actives, 93/35,750 decoys subsample, ratio recorded):
  ROC AUC 0.665, EF@1% 3.33, EF@5% 2.38, BEDROC(α=20) 0.665, with ROC
  plot and per-ligand CSV.

### Scientific honesty (items 4, 5)
- `manifest.receptor.decisions.comparable_to` / `incomparable_to`:
  every run states which preparation engines its scores are comparable
  with; the `none` engine (zero charges/hydrogens) is marked
  incomparable to everything, and the report carries the WARNING.
- `docs/interaction_criteria.md` gains a conventions-vs-heuristics
  table covering every threshold DockFlow applies (success 2.0 Å,
  clustering, score efficiency, contact cutoffs, disulfide/coordination
  detection, box padding and volume guardrails).

### Receptor science (items 6-10)
- `dockflow_core/structure_qc.py`: missing-residue (numbering-gap)
  detection, disulfide detection (S-S < 2.5 Å), reduced-Cys flagging,
  metal-coordination detection (≤ 2.8 Å + geometry labels) - all
  recorded in `manifest.receptor.decisions` and surfaced as warnings.
- `dockflow_core/flexible.py` (ADR-0003): Vina two-file flexible
  side-chain docking - explicit residue lists or contact-derived
  `auto` selection (documented heuristic); score comparability
  degraded and stated.
- Covalent docking via GNINA `--covalent_residue` (ADR-0004); the
  reactive-atom choice stays the researcher's decision.

### Ligand state chemistry (items 11-14, ADR-0008)
- `dockflow_core/ligand_states.py`: protonation (dimorphite-dl, opt-in
  `protonation` extra), tautomers (RDKit), stereocenters, salt
  stripping - default RECORDS the single input state; opt-in
  enumeration docks every state and tags it in the ligand id.

### Analysis (items 16-19)
- Consensus pose ranking (`dockflow rank`, ADR-0005): score z + IFP
  density + geometry term, default weights 0.5/0.3/0.2, always echoed.
- Interaction fingerprints (`dockflow ifp`, ADR-0005): (contact type,
  residue) bitvectors, Tanimoto similarity, greedy leader clustering.

### UX (items 20, 21, 30)
- `visualization/interactive.html`: self-contained 3Dmol.js viewer
  (ADR-0006) with pose switching, grid box, contacts and crystal
  overlay; generated for every completed run (best-effort, never
  fails a run; fixed: viewer mid-run no longer misses the manifest).
- GUI first-run tutorial (item 21): 4-step welcome tour keyed to the
  step bar, QSettings-based first-run detection (version-keyed),
  re-openable from Help; `DOCKFLOW_NO_TUTORIAL=1` for headless runs.
- GUI decision exposure (item 30): the prepare page shows the engine
  comparability sets; the grid-box page shows pocket-provenance
  assumption strength - the same data the manifest records.

### Scale & ops (items 22, 23, 28, 29)
- HPC (ADR-0007): `scripts/generate_sbatch.py` emits a resumable slurm
  array-job template; `dockflow info` detects batch environments and
  warns about cpu oversubscription.
- `scripts/batch_dock.py`: chunked SDF/SMILES streaming (constant
  memory for 100k+ libraries), crash-safe incremental results CSV,
  checkpoint resume, `--backend python` option.
- CI split (item 28): `nightly.yml` runs the full redocking benchmark,
  DUD-E enrichment and the MGLTools comparison on schedule
  (plus `workflow_dispatch`); PR CI stays fast.
- Manifest integrity (item 29): `manifest.checksums` - SHA-256 of every
  key input/output file (run-dir-relative keys); signing deliberately
  not implemented (no key management in an OSS tool - hashes give
  tamper detection).

### Docs
- New ADRs 0003-0008 (flexible residues, covalent docking,
  IFP/consensus, 3Dmol viewer, HPC model, ligand state enumeration).
- Benchmark READMEs rewritten with the real results (the old
  "harness complete, locally validated on 1HVR" status is gone).
- `THIRD_PARTY_LICENSES.md`: 3Dmol.js (BSD-3, CDN) and dimorphite-dl
  (Apache-2.0) recorded; `protonation` extra added to pyproject.
- 1HVR example re-run and committed with full provenance: structure-QC
  decisions, ligand-state provenance, comparability sets, checksums,
  interactive viewer.

## [0.2.1] - 2026-09-08

CI-repair release: fixes the failures observed when the 0.2.0 tree was
pushed to GitHub (test matrix, docker smoke, reproducibility workflow).

### Fixed
- **Docker smoke test exit 127**: the container entrypoint's `shell`
  passthrough forwarded the literal `shell` token to bash
  (`bash shell -c ...` → "No such file or directory"); it now shifts the
  token off first, so `docker run ... shell -c "command -v vina"` works.
- **`broken_openbabel` test fixture crashed on CI matrix cells without
  openbabel** (`ModuleNotFoundError` at setup): the fixture now injects a
  stub `openbabel` package with a valid `__spec__` (so
  `is_importable`/`find_spec` still sees it) when the real wheel is absent
  — macOS, Windows and Py3.10 cells exercise the same broken-dependency
  fall-through as a full install.
- **`test_kabsch_rmsd_matches_bindings_both_paths` hard-required the C++
  accelerator** (unbuilt in the plain test matrix → `ModuleNotFoundError`):
  it now `importorskip`s the module and runs for real in the dedicated
  `bindings-equivalence` CI job, which also executes the dual-path
  kabsch/clustering analyzer tests from now on.
- **`test_pose_cluster_summary_known_structure` asserted translated
  congruent shapes form 3 clusters**: `cluster_poses` measures
  superposition (Kabsch) RMSD, so rigid-body placement cannot separate
  poses — only differing *internal geometry* can.  The test data now uses
  three genuinely distinct conformations (arm-flip variants), and the test
  documents the semantics it verifies.
- **meeko was silently unusable on clean installs**: meeko 0.8.0 imports
  `scipy` (via `meeko.receptor_pdbqt`) without declaring it, exactly like
  its undeclared `gemmi`/`pandas` imports; without scipy every ligand
  preparation fails ("meeko is required for ligand preparation").  `scipy`
  is now pinned in the `prep` and `all` extras — this is also what made
  the `reproducibility` CI workflow (1HVR seed-determinism + scientific
  tests) fail before this release.

### Changed
- `bindings-equivalence` CI job additionally runs the analyzer tests that
  exercise both the C++ and pure-NumPy RMSD paths (`-k "kabsch or cluster"`).
