# ADR 0010: WatVina explicit-water docking backend

- **Status**: proposed (evidence: docs/PORTING_ASSESSMENT.md, Track 1)
- **Date**: 2026-09-08
- **Deciders**: DockFlow maintainers

## Context

No DockFlow stage models displaceable water: "hydration engines" in
the preparator add hydrogens, and Vina scores against a dry receptor.
Water-mediated ligand binding is common enough (bridge H-bonds,
unfavourable displacements) that the review created a
flexible-receptor-and-water track.

WatVina (`biocheming/watvina`) is a self-contained C++ docking engine
(deterministic seed, SDF ligand input, PDB or PDBQT receptors) that
adds: FHFT+GCMC hydration-site prediction (`--predict_water`), water
displacement bias in scoring (`--water_bias_weight`), pharmacophore
generation and template-biased docking, post-docking MM relaxation,
and pybind11 bindings.  License: free for any use ("no license
needed" including non-academic); prebuilt Linux/Windows binaries are
staged in the upstream repo.  It is an *external executable* - same
integration shape as GNINA (ADR-0001), no python dependency.

## Decision

1. Add `docking.backend: "watvina"` following the external-backend
   pattern: executable detection (`DOCKFLOW_WATVINA` env /
   `AppConfig.watvina_exec`), config validation of water options
   (`docking.predict_water`, `docking.water_bias_weight`,
   `docking.water_file`), results parsed from its Vina-style table.
2. Water outputs become first-class artifacts:
   `<run>/analysis/hydration_sites.pdb` (ranked sites, occupancy in
   B-factor), `<run>/docking/water_occupancy.dx` (OpenDX map), both
   linked from report.md.
3. Comparability: WatVina scores (MiniWv model, its own force-field
   terms) are **not comparable** with Vina/GNINA scores; manifest
   gains backend entry + explicit `incomparable_to` list.  Receptor
   input flavour (raw PDB vs prepared PDBQT) is recorded because
   charge models then differ.
4. Not installed by any extra; `scripts/install_tools.sh` gains an
   opt-in downloader for the upstream release binary, gated on
   confirming redistribution terms with the author (action item in
   PORTING_ASSESSMENT.md); CI uses fake-executable tests only.
5. Validation before promotion out of experimental: rerun the
   24-complex redocking benchmark with water prediction on; publish
   pose-recovery deltas with and without water bias.

## Consequences

- The water track becomes testable science rather than a roadmap
  bullet: benchmark evidence or it stays experimental.
- One more backend in the manifest/env fingerprint surface - the
  comparability matrix grows a row; users must not mix backends
  mid-screening (documented).
- Risk: single-maintainer upstream with a moving CLI surface - pin
  the tested engine version in `AppConfig` defaults and record the
  actual binary's `--version` per run.
