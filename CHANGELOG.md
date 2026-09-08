# Changelog

## [1.0.0rc1] - 2026-09-08

**Feature freeze for v1.0 release and paper submission.**  Final
scientific-rigor tightening pass (20 review items): no new features, no
new pipeline stages, no new runtime dependencies - every change either
makes a published claim more honest, adds a missing statistic, or makes
a published number regenerable.  The post-freeze roadmap lives in
`ROADMAP_POST_V1.md`.

### Fixed - CI (all three gui jobs red on the rc1 push)
- **The `qapp` fixture was module-private.**  The new final-pass test
  `test_confirm_degraded_engine_persists_acknowledgement` (item 8, GUI
  modal persistence) requested the `qapp` fixture, which was defined
  inside `test_gui_smoke.py` only - invisible to every other module, so
  the test errored at setup ("fixture 'qapp' not found") and took the
  gui job down on all three operating systems.  The fixture now lives
  in `tests/conftest.py` (session-scoped; imports PyQt6 lazily so
  non-GUI runs never need it; skips gracefully when the binding or its
  system GL runtime is unavailable), and `test_gui_smoke.py` consumes
  the shared fixture instead of defining its own.
- **Two analyzer tests failed in RDKit-less environments.**
  `test_crystal_rmsd_with_method_prefers_rdkit_path` and
  `test_crystal_rmsd_getbestrms_benzoate_symmetry` assert that the
  RDKit `GetBestRMS` graph-automorphism path is chosen; in minimal
  environments (the gui CI job installs `[test,gui,viz]` without the
  `prep` extra) RDKit is absent and the analyzer legitimately falls
  back to element-greedy Kabsch, so the assertions failed instead of
  skipping.  Both tests now carry an explicit `requires_rdkit` skipif
  guard, matching the established pattern in test_batch_dock and
  test_preparator.

### Fixed - scientific correctness (found BY the new validation)
- **Symmetry-aware crystal RMSD was atom-order dependent.**  The
  analyzer seeded the pose-to-reference correspondence geometrically and
  then let `GetBestRMS` permute only over graph automorphisms - when the
  pose file and the reference PDB listed atoms in different orders (which
  Vina/meeko outputs and crystal PDBs routinely do), the seed could lock
  in a topologically wrong assignment that the automorphism search could
  not repair, over-penalising symmetric ligands by up to ~0.5 Å.  Found
  by the new independent spyrmsd cross-check (item 5); fixed by building
  the pose molecule with its own graph so `GetBestRMS` minimises over
  all graph isomorphisms.  After the fix: **193 benchmark poses,
  max |Δ| vs spyrmsd = 0.0 Å**.  The same bit-identical poses were
  re-analysed with the corrected RMSD, so the published benchmark
  numbers changed (they were conservative, not optimistic): best-any-rank
  success 18/24 → **23/24**, rank-1 (the new formal definition)
  **21/24**, mean best-pose RMSD 1.32 Å → **0.67 Å**; five of the six
  v0.3.x "failures" actually recover a ≤ 2.0 Å pose.  The remaining
  failure modes are now classified (ranking vs sampling) via the
  published top-N curve.
- comparison.csv silently dropped `mgl_rmsd_a` (and three other
  MGL-arm columns) from the published file although the summary table
  used them - the full column set is now written.
- `dockflow run --workdir` was accepted but ignored (the YAML workdir
  always won).
- single-stage resume re-saved manifests WITHOUT the provenance sections
  (receptor decisions, grid box, ligand prep, docking backend) of the
  stages that did not re-run; they are now restored and merged.
- report.md: double period after the RMSD method note; stale
  hard-coded fallback version string in source checkouts.

### Changed - claim honesty (the friend's final review)
- **Pose recovery has a formal definition** (rank-1 pose, symmetry-aware
  heavy-atom RMSD ≤ 2.0 Å, community convention per Trott & Olson 2010 /
  CSAR / D3R practice, NOT a validated cutoff) and the benchmark
  publishes ALL definitions (rank-1, top-3, legacy best-any-rank) at
  1.0/2.0/3.0 Å plus the top-N success curve - the number cannot be
  definition-shopped (docs/interaction_criteria.md).
- **Protonation honesty**: "pH 7.4 convention" is gone - receptors keep
  DEPOSITED states + template hydrogens (pH 7.4 nominal, NOT a pKa
  calculation); the report footer says so, and a new catalytic-residue
  warning names the His/Asp/Glu residues near the box with PROPKA/H++
  guidance.  README + docs/preparation_assumptions.md (new, the full
  does/does-not/should table for every preparation assumption) rewritten
  accordingly.
- **Degraded mode is unmistakable** (item 8): the `none` engine now
  triggers a ⚠ SCIENTIFICALLY DEGRADED banner at the top of report.md,
  a blocking GUI modal ([Cancel run] / [I understand, continue anyway],
  persisted per version in QSettings), a red-stderr CLI gate requiring
  `--allow-degraded` to proceed, and machine-readable
  `manifest.receptor.decisions.degraded_mode/degraded_reason`.
- **Grid-box provenance in every report** (item 9): the box line names
  the exact pocket atom (`pocket:XK2 chain A residue 263`), padding and
  `assumption: strong|medium|weak|user-explicit`; weak/medium boxes and
  >60%-coverage boxes carry warnings (wired since 0.3.0, now surfaced in
  the report line as well).
- **README claim audit** (item 16): beta banner ("no claim of superior
  performance"), "MGLTools equivalent" → "MGLTools-equivalent *logic*
  (NOT the same code, NOT bit-identical)", "less accurate fallback" →
  "scientifically degraded", "validated example" → "evaluated
  single-complex example", every "automatic" qualified as mechanical;
  a CI doc-lint fails on unqualified flagged words, and
  docs/claim_audit.md records every change (item 17).
- **Enrichment reframed as an evaluation** (item 11): the AUC line now
  states single-target/exhaustiveness-1/no-comparative-baseline
  explicitly; the BEDROC random baseline (NOT 0.5 at α=20) is measured
  over 10 000 seeded permutations instead of assumed.

### Added - benchmark completeness and statistics
- Per-complex preparation + grid provenance columns in
  `redocking_results.csv` (item 1: receptor engine/atoms, ligand
  rotatable bonds/heavy atoms, grid box source/center/size/padding/
  assumption strength, exhaustiveness/modes/seed/cpu, RMSD method) -
  a reader can now tell whether a failure is prep- or docking-related.
- Wilson score 95% CIs for every success rate (hand-rolled, formula
  documented; 18/24 example verified as 55.1%-88.0%), RMSD distribution
  (min/p25/median/p75/max/std) + histogram PNG, per-target grouping
  sub-table, failed-complex sub-tables with hypothesized causes and
  failure-mode classification (item 2).
- **MGLTools comparison fully documented** (item 4): full crash
  tracebacks in `results/failures/*.txt` (the 200-char CSV truncation
  is gone), `failure_stage` + `failure_root_cause` columns, the explicit
  shared-denominator statement, Wilson CIs for both arms; both arms
  re-run with the corrected RMSD (DF 20/21 vs MGL 18/21).
- **DUD-E enrichment on 3 targets** (item 12): aa2ar (3EML) and parp1
  (2RD6) added with the identical seeded protocol; committed
  SHA-256-pinned ISM lists; `--reproduce` re-docks the exact 133-ligand
  panel and asserts the metrics within 1e-6 (item 10); duplicate-SMILES
  verification documented (7 duplicates in the hivpr panel, reported
  not dropped); scoring/tie-breaking/metric formulas documented in the
  benchmark README.
- **Regeneration harness** (item 14): `complexes.csv` pins resolution +
  UniProt + SHA-256 of every committed raw PDB;
  `benchmarks/redocking/regenerate.py` re-runs the full benchmark with
  the thread count pinned (Vina is thread-count dependent - the v0.3.1
  CI lesson) and ASSERTS every published number within 1e-6 (verified:
  0 mismatches); the nightly workflow runs it and fails on drift.
- **Symmetry-RMSD validation** (item 5): analytic known-answer tests
  (benzene 60°, naphthalene flip, aspirin O-swap, chiral ligand,
  1.5 Å-displaced pose) + the spyrmsd cross-check script with
  `rmsd_crosscheck.png` and a <0.1 Å acceptance gate; the
  `crystal_rmsd_method` column in the benchmark CSV.
- **Preparation validation extended** (item 13): 24-complex
  DockFlow-vs-MGLTools aggregate + meeko-CLI faithfulness check,
  generated from committed artifacts by
  `benchmarks/preparation_validation_aggregate.py`.
- **Paper figures** (item 19): `paper/figures/fig1-5.{png,pdf}` rendered
  from the committed benchmark CSVs only, with a provenance README.
- **docs/SCIENTIFIC_VALIDATION.md** (item 20): the single citable page
  (7 sections: benchmarks, comparison, enrichment, RMSD validation,
  preparation validation, limitations, regeneration).

## [0.3.1] - 2026-09-08

CI-stability release: fixes both failing GitHub Actions jobs
(windows tests + scientific tests) and adds the upstream-reuse
assessment requested for roadmap planning.

### Fixed - CI (the two red jobs)
- **windows-latest tests (1 failure)**: the sbatch bash-syntax test
  invoked `bash` on PATH, which on Windows resolves to the WSL
  launcher stub that prints a UTF-16 "no distro" error and exits 1.
  `scripts/generate_sbatch.py` now writes LF-only templates
  (`newline="\n"` - sbatch targets Linux clusters; CRLF would break
  bash/slurm anyway), the test writes bytes + asserts LF endings +
  shebang on every OS, and `bash -n` runs only where a real POSIX
  bash exists (Git-for-Windows bash on Windows, skip otherwise).
  New always-on test `test_template_is_lf_only`.
- **scientific tests (RMSD 2.54 A > 2.0 A on a 4-vCPU runner)**: the
  1HVR pose-recovery assertion was machine-dependent.  Vina results
  depend on thread count (same seed, different cpu count -> different
  search trajectory), and `receptor.engine: auto` resolved to rdkit
  in CI (no openbabel there) but openbabel locally - two silent
  divergences.  The test config is now fully pinned: `engine=rdkit`
  (available on every platform via the prep extra), `cpu=1`
  (deterministic on any core count), exhaustiveness 2 -> 8 (default
  search effort).  Validated locally: best pose crystal RMSD 0.54 A,
  best affinity -11.20 kcal/mol; the misleading "deterministic
  regardless of thread count" docstring is corrected.  The CI job
  pins vina/meeko/rdkit versions (validation job = known-good env).
- Corrupted branch filters (`branches: ain]` -> `branches: [main]`)
  fixed in build.yml and reproducibility.yml (the corruption had
  been fixed before in 0.2.0 and regressed in the v0.3.0 upload).

### Added - upstream-reuse assessment (ADR evidence)
- `docs/PORTING_ASSESSMENT.md`: twelve proposed upstream projects
  (WatVina, RxDock, OpenPharmacophore, PharmacoNet, OpenDock,
  PandaDock, ProLIF, PLIP, 3Dmol.js, NGLview, Molecular Nodes,
  Avogadro) reviewed for portable functions/features with verified
  licenses, integration routes (dep / external binary / vendored /
  ideas-only) and a recommended sequence.
- `ADR-0009` (proposed): ProLIF optional IFP backend - richer
  interaction vocabulary + per-pose LigNetwork 2D diagrams (the new
  Track 2), homegrown IFP stays the zero-dependency default.
- `ADR-0010` (proposed): WatVina external-binary backend - explicit
  water (FHFT+GCMC hydration sites, displacement bias) docking with
  comparability caveats, same pattern as the GNINA backend.

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
