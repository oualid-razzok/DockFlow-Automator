# ADR 0009: ProLIF as optional interaction-fingerprint backend

- **Status**: proposed (evidence: docs/PORTING_ASSESSMENT.md, Track 2)
- **Date**: 2026-09-08
- **Deciders**: DockFlow maintainers

## Context

DockFlow computes interaction fingerprints (IFPs) with a homegrown,
deliberately small vocabulary (analyzer `ifp_bits` /
`ifp_vocabulary`) used by the consensus ranker (ADR-0005).  The
second-round review asked for 2D ligand-interaction diagrams as a new
roadmap track; a diagram needs richer interaction semantics
(donor/acceptor-direction H-bonds, pi-stacking, salt bridges, vdW)
than our geometric heuristics provide.

ProLIF (Apache-2.0, actively maintained, on PyPI/conda-forge) is the
community-standard library for exactly this: RDKit-based interaction
fingerprints for protein-ligand complexes from docking poses, plus
**LigNetwork** browser-rendered 2D interaction diagrams.  The standing
dependency constraint (v0.2.1 work order) allows new optional
dependencies only with an ADR - this is that ADR.

## Decision

1. Add an optional `prolif` extra (`prolif>=2.0`); never a hard or
   `all` dependency - the homegrown IFP remains the zero-dependency
   default.
2. When `prolif` imports, analyzer computes the canonical ProLIF
   fingerprint for every analysed pose (per interacting residue pair)
   and writes a per-pose 2L LigNetwork diagram to
   `<run>/visualization/lignetwork_<ligand>_<pose>.html` - same
   self-contained-HTML philosophy as the 3Dmol viewer (ADR-0006).
3. Fingerprint provenance: the manifest records
   `analysis.ifp_flavour: "homegrown" | "prolif"`; consensus ranks
   produced from different flavours are flagged as not comparable
   (extends the `comparable_to` machinery, item 4).
4. The consensus IFP term (ranking.py) uses the ProLIF fingerprint
   when present; weights unchanged and still reported per component.
5. ProLIF deps (rdkit only for static structures) overlap our `prep`
   extra - no new transitive burden beyond `prolif` itself.

## Consequences

- Users get publication-grade 2D diagrams without DockFlow adopting a
  rendering stack: ProLIF emits the D3 page, we frame and link it.
- Two fingerprint flavours to keep honest: tests must cover the
  adapter (monkeypatched prolif) plus a flavour-provenance test;
  vocabulary drift in upstream ProLIF is a documented risk recorded
  in the environment fingerprint.
- Apache-2.0 attribution added to THIRD_PARTY_LICENSES.md (NOTICE
  style), and the license-audit script learns the new optional dep.
- If ProLIF is absent: zero behaviour change (skips silently, one
  INFO log line, report notes "2D diagrams need the prolif extra").
