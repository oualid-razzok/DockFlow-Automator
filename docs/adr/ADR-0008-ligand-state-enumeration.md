# ADR 0008: Ligand state enumeration (protonation, tautomers, stereo, salts)

- **Status**: accepted (work-order items 11-14)
- **Date**: 2026-09-08

## Context

A docking input SMILES/SDF fixes ONE protonation state, ONE tautomer,
ONE stereochemistry - silently.  Different states can dock completely
differently.  The work order ordered the chemistry to be surfaced:
detected, reported, optionally enumerated, never silently changed.

## Decision

1. **Default = record, do not change**: the manifest's
   `ligands[*].prep` block records the input's protonation source
   (`dimorphite-dl` pH calculation or "input as-given"), tautomer
   choice, undefined stereocenters and salt fragments - so the
   researcher knows which ONE state was docked.
2. **Opt-in enumeration**: `protonation: enumerate` docks every
   microspecies dimorphite-dl considers populated at the target pH
   (7.4 +/- 0.5 default); `tautomers: enumerate` uses RDKit's
   `TautomerEnumerator` (canonical set, NO pH weighting - limitation
   recorded); `stereo: enumerate` docks every stereoisomer of
   undefined centers; salts are stripped with an exact fragment
   report (what was removed, what was kept).
3. Each enumerated state becomes its own PDBQT and its own
   `summary.csv` row (ligand id carries the state tag), so no state is
   ever silently merged into another.
4. `dimorphite-dl` is an optional dependency (allowed by the work
   order without further process); without it, protonation
   enumeration is an actionable error, and pH-based reporting degrades
   to "input as-given".

## Consequences

- Enumeration multiplies docking cost by the number of states - it is
  opt-in, and the manifest records how many states each ligand
  produced.
- RDKit tautomer enumeration is rule-based (not population-weighted);
  the manifest says so, and results are labelled accordingly.
