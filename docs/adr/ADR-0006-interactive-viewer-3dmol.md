# ADR 0006: Interactive 3D viewer via 3Dmol.js

- **Status**: accepted (work-order item 20)
- **Date**: 2026-09-08

## Context

PyMOL renderings are static images; researchers want to rotate/zoom
poses and inspect contacts interactively.  A shipped GUI 3D engine
(PyQt+OpenGL widget) would add heavy dependencies and platform
fragility for marginal benefit.

## Decision

1. Every completed run gets a **self-contained HTML viewer**
   (`visualization/interactive.html`) built on **3Dmol.js** loaded from
   its CDN.  This is the one new front-end dependency explicitly
   allowed by the work order (no build step, no npm lockfile).
2. The HTML embeds all data at generation time (receptor + poses as
   inline PDBQT text, grid box, contacts, crystal reference overlay
   when redocking validation exists) - the file works offline after
   generation except for the CDN script itself.
3. Pose switching, box display, contact highlighting and the crystal
   overlay are client-side JavaScript; no server component, no
   tracking, no external data flow beyond the CDN fetch.
4. The viewer is generated at the end of the visualization stage
   (after a manifest snapshot is published) and is best-effort: a
   viewer failure NEVER fails the run.

## Consequences

- Opening the viewer requires internet (CDN) but no installation;
  fully air-gapped users keep the PyMOL/matplotlib static renders.
- `THIRD_PARTY_LICENSES.md` records 3Dmol.js (BSD-3) and the CDN
  origin.
