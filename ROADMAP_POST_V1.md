# Roadmap after v1.0 (feature-frozen scope)

**Status: v1.0.0-rc1 is FEATURE-FROZEN for the v1.0 release and the
paper submission.**  Nothing below this line is in the paper, in the
v1.0 milestone, or needed to reproduce any published number.  Each entry
records WHY it was deferred (not merely that it was), and the porting
assessment for upstream ideas lives in `docs/PORTING_ASSESSMENT.md`
(license-verified, per-repo routes; ADR-0009 and ADR-0010 exist for
ProLIF and WatVina).

## Deferred scientific features

| feature | why deferred | entry point |
|---|---|---|
| ProLIF interaction fingerprints (Track 2 LigNetwork 2D diagrams) | the v1.0 IFP (bit-per-residue fingerprints) is committed and benchmarked; the richer ProLIF backend is an optional dependency and a new validation surface | ADR-0009 (Apache-2.0, optional extra) |
| WatVina explicit-water docking | a new backend + water-aware validation benchmark does not exist yet; must not enter v1.0 unbenchmarked | ADR-0010 (external binary) |
| RxDock port (flexible-receptor algorithms) | LGPL-2.1 - process-boundary port only; large engineering surface with no validation budget before the paper | `docs/PORTING_ASSESSMENT.md` |
| AD4 maps (`autogrid`) as a first-class scoring path | the partial-charge differences documented in `docs/preparation_validation.md` must be handled with a dedicated engine-comparability note first | — |
| Ensemble / multi-receptor docking + consensus rescoring | implemented primitives exist (consensus rank column), but no ensemble benchmark; publishing unbenchmarked features is what the final pass removed | — |
| 2D ligand-interaction diagrams | Track 2 (post-paper); PLIP is GPL-2.0 → CLI-boundary only | `docs/PORTING_ASSESSMENT.md` |
| Plugin API for custom scoring (Smina custom scores) | API stability commitment belongs after the paper freeze | — |
| Covalent / metal-aware scoring benchmarks | covalent docking is implemented but unbenchmarked - it stays listed in the limitations of `docs/SCIENTIFIC_VALIDATION.md` until a benchmark exists | — |

## Deferred validation work

| study | why deferred |
|---|---|
| Full DUD-E (102 targets) virtual screening | compute budget; the 3-target evaluation states its own smallness |
| Multi-seed variability study | one seed everywhere is a stated limitation; measuring the spread is a paper of its own |
| Cross-docking benchmark (ligand into a different structure of the same target) | needs curated target panels; the 1HSG example config exists but no published rate |
| Flexible-receptor benchmark | flexible residues are implemented and flagged incomparable; no benchmark yet |
| Higher-exhaustiveness / scoring-function comparison grid | exhaustiveness 8 is pinned for the published numbers; a grid changes the sampling-failure boundary and must be rerun as a whole |

## Not planned (explicit non-goals)

* **QM-derived charges** (RESP/AM1-BESP) - preparation stays
  Gasteiger-with-provenance; QM rescoring belongs to downstream tools.
* **pKa calculation inside the pipeline** - PROPKA/H++ remain external
  inputs by design: the pipeline records that it does NOT do pKa
  (see `docs/preparation_assumptions.md`) rather than approximating it.
* **New runtime dependencies before v1.1** - dimorphite-dl (protonation
  enumeration) and NGL/3Dmol.js (viewer) are the only optional extras
  added since 0.2.x; anything else requires an ADR.
