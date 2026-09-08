# ADR 0005: Interaction fingerprints (IFP) and consensus pose ranking

- **Status**: accepted (work-order items 18 and 19)
- **Date**: 2026-09-08

## Context

Ranking poses by Vina affinity alone ignores binding-mode plausibility
(the classic "best score, wrong pose" failure).  Two work-order items
requested (a) ligand-receptor interaction fingerprints with pose
similarity/clustering, and (b) a consensus pose ranking that combines
score, interaction pattern and geometry.

## Decision

1. **IFP** (`dockflow ifp`, `dockflow_core/analyzer.py`): each pose is
   encoded as a bit vector over (contact type, residue) pairs derived
   from the EXISTING geometric contact analysis - no new criteria are
   invented (same cutoffs as `docs/interaction_criteria.md`).  Pose
   similarity = Tanimoto; clustering = greedy leader clustering at a
   configurable similarity threshold (default 0.6).  RDKit
   `ExplicitBitVect` when available, pure-python fallback otherwise.
2. **Consensus ranking** (`dockflow rank`, `dockflow_core/ranking.py`):
   `score = w_score * z(score) + w_ifp * IFP_density + w_rmsd *
   geometry_term` with default weights 0.5 / 0.3 / 0.2 - all three
   documented as **tunable heuristics with no thermodynamic meaning**.
   The rank output labels itself accordingly and prints the weights.
   When IFP or geometry data are missing, the ranking degrades to
   score-only and SAYS so (never silently).
3. Both operate on COMPLETED runs (manifest.json + docking outputs) -
   they are post-hoc analysis commands, not pipeline stages, so the
   core pipeline stays deterministic and replayable.

## Consequences

- No new mandatory dependencies; bit-vector performance is fine for
  the 9-20 poses per ligand that DockFlow produces.
- Consensus ranks change when weights change; every output therefore
  echoes the weights used, so a published table can be reproduced.
