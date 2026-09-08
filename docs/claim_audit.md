# Claim audit (final pass, items 16-17)

Every occurrence of the five review-flagged words — **equivalent,
validated, accurate, automatic, production-ready** — across the v1.0.0-rc1
documentation, with the qualification applied and the evidence behind it.
This file is itself auditable: the CI doc-lint
(`tests/test_final_pass.py::test_flagged_claims_are_qualified`) fails if
any flagged word appears unqualified in the linted docs, and
`test_claim_audit_document_exists` requires this file to record the
decisions.

## "equivalent" (MGLTools)

| Where | Original claim | Qualification applied | Evidence |
|---|---|---|---|
| README pipeline-stages table | "MGLTools equivalent" | "MGLTools-equivalent *logic* (reimplemented on Meeko+RDKit+OpenBabel; NOT the same code, NOT bit-identical output)" | `docs/preparation_validation.md` (per-property comparison + 24-complex aggregate: receptor heavy-atom counts identical, torsion perception differs on 6/21, charge sums differ on ionised groups) |
| README repository layout | "MGLTools-equivalent receptor/ligand preparation" | "MGLTools-*logic* … (reimplemented …; NOT bit-identical)" | same as above |
| CHANGELOG 0.2.0-era wording | historical entries kept as written at the time; current entries use the qualified form | — | git history |

## "validated"

| Where | Original claim | Qualification applied | Evidence |
|---|---|---|---|
| README "Validated example: 1HVR redocking" | implied a validated pipeline | "Validated example … (single complex, 1 of 24)" — validated in the narrow sense of *one complex evaluated*; full suite in benchmarks/ | `examples/results/1HVR/report.md` (RMSD 0.84 Å), `benchmarks/redocking/results/summary.md` (24 complexes, all definitions + CIs) |
| README repository layout | "artifacts of the validated 1HVR run" | "artifacts of the 1HVR example run (single complex; the 24-complex suite is in benchmarks/)" | same |
| docs/interaction_criteria.md | "validated thresholds" table heading | kept, but the table itself classifies every threshold as Convention/Heuristic, explicitly "not validated for any specific target" | the table itself |

## "accurate"

| Where | Original claim | Qualification applied | Evidence |
|---|---|---|---|
| README degradation note | "a working (if less accurate) fallback, never a crash" | "a working — but **scientifically degraded** (geometry-only, zero charges; CLI requires `--allow-degraded`, GUI asks, report carries the ⚠ banner) — fallback, never a silent crash" | `dockflow_core/preparator.py` (none engine), `degraded_mode_banner` in pipeline.py, tests/test_final_pass.py |

## "automatic"

The principle (recorded in the README): *automatic engine selection is a
mechanical choice; the scientific choices — protonation, tautomer,
stereochemistry — are recorded in the manifest and are the researcher's
to override.*

| Where | Original claim | Qualification applied | Evidence |
|---|---|---|---|
| README features (backend selection) | "selected automatically" | "selected automatically (a mechanical availability choice; the scientific choices … stay yours and are recorded)" | `dockflow_core/docker_engine.py` + manifest |
| README features (matplotlib fallback) | "automatic matplotlib fallback" | "automatic (mechanical renderer choice, not a scientific one)" | `dockflow_core/visualizer.py` |
| README quick start (vina CLI install) | "installed automatically by the scripts" | "installed automatically by the scripts below — a mechanical install step, no scientific choices involved" | `scripts/install_tools.*` |
| README quick start (engine download) | "automatic Vina engine download" | "automatic (mechanical) Vina engine download" | same |
| README provenance section | '"automatic preparation" must never mean "hidden preparation"' | kept — this line IS the qualification | `manifest.json` decisions blocks |
| USER_GUIDE / BUILD_GUIDE | "automatically" used for install/caching mechanics | reviewed; all remaining uses describe mechanical steps (downloads, caching, retries), none describe scientific decisions | docs |

## "production-ready"

| Where | Original claim | Qualification applied | Evidence |
|---|---|---|---|
| (none found in the linted docs) | — | the phrase is banned outright; the degraded-mode banner instead says "NOT recommended for production" | `degraded_mode_banner` in pipeline.py |

## Other honesty changes made in the same pass (not word-flagged, recorded for completeness)

* README beta banner added (BETA status, validation-suite pointer, "no
  claim of superior performance").
* The default-protonation bullet rewritten (deposited states + template
  hydrogens at pH 7.4 nominal, NOT a pKa calculation, PROPKA/H++
  guidance) with cross-links from the report footer —
  `docs/preparation_assumptions.md` carries the full table.
* Benchmark headline numbers re-qualified to the formal rank-1 pose
  recovery definition with Wilson 95% CIs, plus the alternative
  definitions (top-3, best-any-rank) —
  `benchmarks/redocking/results/summary.md`.
* Enrichment AUC reframed as an EVALUATION (not a superiority claim) with
  random baselines measured explicitly —
  `benchmarks/enrichment/results/<target>/summary.md`.
* The v0.3.x "18/24 (75%)" success rate was corrected to
  "rank-1 21/24, top-3 23/24, best-any 23/24" after the symmetry-RMSD
  atom-order fix (see CHANGELOG 1.0.0rc1: the old numbers were
  over-penalised, not over-claimed).
