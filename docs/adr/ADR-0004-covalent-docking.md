# ADR 0004: Covalent docking via GNINA

- **Status**: accepted (work-order item 9)
- **Date**: 2026-09-08

## Context

Covalent inhibitors (irreversible binders, e.g. cysteine-reactive
fragments) cannot be redocked meaningfully with standard Vina: the
warhead-ligand bond does not exist in the search.  The GNINA backend
supports covalent docking natively
(`--covalent_residue CHAIN:RESNUM:ATOM` + `--covalent_lig_atom_idx` +
`--covalent_box_padding`, default reactive bond length 1.7 A).

## Decision

1. Covalent docking is GNINA-only.  `docking.covalent_residue` /
   `docking.covalent_ligand_atom` config keys map onto the GNINA
   flags; requesting them with the vina backend is a config error with
   a redirect to GNINA (never silently ignored).
2. The box for a covalent run is auto-centered on the target residue
   when the researcher does not supply one (weak assumption, recorded
   as such in `manifest.gridbox.assumption_strength`).
3. GNINA is optional (`app.gnina_exec`); its absence produces an
   actionable error, not a crash.  No attempt is made to emulate
   covalent docking with Vina (dishonest).
4. Provenance: `manifest.docking.covalent` records residue, ligand
   atom index and bond length; the report labels the run "covalent
   (GNINA)" and warns that CNN scoring and covalent scoring are NOT
   comparable to standard Vina scoring.

## Consequences

- Covalent workflows require a GNINA installation (not bundled);
  everything else (receptor/ligand prep, analysis, RMSD) is shared
  with the non-covalent path.
- The reactive-atom choice is the researcher's scientific decision -
  DockFlow records it, never chooses it.
