# ADR 0002: Docked-pose SDF export

- **Status**: accepted (work-order item 29)
- **Date**: 2026-09-07

## Context

Vina writes multi-pose PDBQT; most chemistry tooling cannot read PDBQT.
The work order requested `<ligand>_poses.sdf` with one record per pose
and SD fields for the scores and crystal RMSD.  This introduces a new
*output* file format (the work order itself requested it; this ADR
records the design).

## Decision

1. `docking.export_sdf: true` (default) in the pipeline config; can be
   disabled.  Not a new pipeline stage: the export runs inside the
   analysis stage (after `crystal_rmsd` is known, so SD fields can
   include it).
2. Conversion strategy, most rigorous first:
   Meeko `PDBQTMolecule` + `RDKitMolCreate` (restores bond orders from
   PDBQT remarks) → per-pose RDKit `MolFromPDBBlock` with proximity
   bonding as fallback.
3. SD fields: `vina_affinity`, `vina_rmsd_lb`, `vina_rmsd_ub`,
   `crystal_rmsd` (when a redocking reference exists), plus
   `docking_backend`.
4. Failures are warnings, never fatal (SDF export is a convenience
   layer on top of the authoritative PDBQT).

## Consequences

- Outputs are meeko-dependent but degrade gracefully without it.
- The docking directory now contains PDBQT (authoritative), SDF
  (interoperability) and CSV/JSON (analysis) representations of the
  same poses.
