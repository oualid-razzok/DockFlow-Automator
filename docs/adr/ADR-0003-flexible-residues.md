# ADR 0003: Flexible receptor residues (Vina two-file side-chain docking)

- **Status**: accepted (work-order item 8)
- **Date**: 2026-09-08

## Context

Rigid-receptor redocking fails exactly on the complexes where the
binding site adapts to the ligand (the benchmark's six failures are
all flexible peptidomimetic HIV-PR inhibitors).  AutoDock Vina supports
side-chain flexibility through a two-file receptor format: a rigid
PDBQT plus a flex PDBQT with `BEGIN_RES`/`END_RES` torsion-tree blocks.
The legacy path was MGLTools' `prepare_flexreceptor4.py`.

## Decision

1. `dockflow_core/flexible.py` generates both files from an ALREADY
   prepared rigid receptor PDBQT plus a residue selection
   (`receptor.flexible_residues: ["A:ASP25", ...]` or `auto`).
2. Rotatable side-chain bonds come from a per-residue-type chi table
   (chi1..chi4 of the standard side-chain flexibility definitions) -
   deterministic, matches `prepare_flexreceptor4.py` for the common
   cases, no ring/aromaticity perception pass required.
3. The flexible part is the side chain from CB (SG for CYS) onward;
   backbone (N, CA, C, O, OXT) always stays rigid.  GLY/ALA/PRO have no
   rotatable side-chain bonds - requesting them warns, never errors.
4. `auto` selection derives the residue list from geometric contacts
   between the co-crystallized reference ligand and the receptor
   (residues with any contact become flexible) - a documented,
   overridable heuristic, not a simulation.
5. **Score comparability is degraded and the run says so**: flexible-
   receptor scores are NOT comparable with rigid-receptor scores of the
   same ligand (the flexible side chain pays an entropy-like penalty
   inside the Vina scoring); the manifest and report carry an explicit
   warning and `manifest.receptor.decisions.flexible_residues`.

## Consequences

- Hard redocking cases become addressable at the cost of a heuristic
  chi-table (non-standard residues fall back to "largest rotatable
  side-chain subtree" with a warning).
- The two files must be regenerated whenever the residue selection or
  the rigid preparation changes; `manifest.paths` records both.
