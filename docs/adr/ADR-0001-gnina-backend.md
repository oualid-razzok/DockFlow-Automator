# ADR 0001: GNINA docking backend

- **Status**: accepted (work-order item 32, explicitly permitted new
  optional dependency)
- **Date**: 2026-09-07
- **Deciders**: DockFlow maintainers

## Context

GNINA (Apache-2.0) is the community-standard ML-augmented docking/scoring
engine.  Users asked for CNN rescoring and CNN-guided docking alongside
Vina/Smina.  GNINA is an *external executable*, not a python package, so
adding it introduces no python dependency and no redistribution concern.

## Decision

1. Add `docking.backend: gnina` (CLI auto-fallback chain:
   python-vina → vina CLI → smina → gnina).
2. GNINA-specific options, validated at config load:
   `docking.cnn_scoring: rescore|all|none` and `docking.cnn: <model>`
   (rejected unless `backend: gnina`).
3. GNINA reuses the Vina-style long options it shares (receptor, ligand,
   out, center/size, exhaustiveness, num_modes, seed, cpu); its log goes
   to stdout (no `--log`); the Vina results table it prints is parsed by
   the existing `parse_vina_log`.
4. `DOCKFLOW_GNINA` env var / `AppConfig.gnina_exec` for explicit paths,
   mirroring vina/smina.
5. Not installed by any extra; documented install hint
   (conda-forge gnina, or the upstream release binaries).

## Consequences

- `detect_backends()` and `dockflow info` report GNINA availability.
- Scoring semantics differ from Vina (CNN): absolute scores are not
  comparable across backends; the manifest records the resolved backend
  per run (stage audit trail + environment fingerprint).
- Tests use a fake cross-platform executable (same pattern as the fake
  vina CLI tests); no real GNINA run in CI (binary size/licensing of
  CNN weights is the user's environment decision).
