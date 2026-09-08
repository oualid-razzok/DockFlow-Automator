# Upstream reuse & porting assessment

Reviewed 2026-09-08 against v0.3.1 planning.  Scope: the twelve upstream
projects proposed as sources of portable functions and features, grouped
by track.  Each entry answers four questions: **what it has that DockFlow
lacks, what exactly we would port, what the license allows, and what the
integration route is** (optional dependency / external binary backend /
vendored code / ideas-only).

DockFlow is MIT-licensed.  The standing dependency constraint (v0.2.1 work
order) is: *no new pipeline stages, file formats or optional dependencies
beyond `dimorphite-dl` and NGL/3Dmol.js unless an ADR is written*.  This
document is the evidence base for those ADRs; it proposes two (ADR-0009,
ADR-0010) and explicitly rejects or defers the rest with reasons.

## License compatibility rules applied

| Upstream license | Import as Python dep | Vendor code into repo | Use as external binary/CLI |
|---|---|---|---|
| MIT, BSD-3, Apache-2.0 | yes (attribution) | yes (attribution + notice) | yes |
| LGPL-2.1 | no (import implies LGPL obligations) | no | yes (process boundary) |
| GPL-2.0 | **never** | **never** | yes (process boundary, like PyMOL) |
| custom/academic (WatVina) | n/a (C++ binary) | no | yes, verify redistribution terms |

The process-boundary rule is not new: DockFlow already runs Vina, GNINA
and PyMOL as external executables, and ADR-0001 (GNINA) established the
pattern - external engines are detected, named in the manifest, and their
results carry comparability caveats.  GPL/LGPL code may never be linked
or imported, only executed.

## Summary table

| Repo | License | Route | Priority | One-line verdict |
|---|---|---|---|---|
| WatVina | custom (free use) | external binary backend | **high** | explicit-water docking + pharmacophore bias in one self-contained binary - the missing water track |
| ProLIF | Apache-2.0 | optional Python dep | **high** | richer IFP chemistry + 2D LigNetwork diagrams (Track 2) behind an adapter |
| PLIP | GPL-2.0 | external CLI/Docker only | medium | best interaction *reporting*; never import - GPL is viral |
| OpenDock | MIT | port patterns, not code | medium | constraint objects + scoring-function registry are the portables |
| PandaDock | MIT | optional Python dep (engine) | medium | pure-Python engine = Windows/macOS fallback + GNN rescoring for consensus |
| RxDock | LGPL-2.1 | external binary backend | low/defer | flex receptor + pharmacophore restraints, but WatVina covers more with less build pain |
| PharmacoNet | MIT | ideas only | low | DL pharmacophores are excellent but torch+pymol deps are too heavy |
| OpenPharmacophore | MIT | vendor if track starts | low | feature-perception code is clean; project dormant since 2023 |
| 3Dmol.js / py3Dmol | BSD-3 | already used; vendor JS + notebook embed | **high** (small) | local 3Dmol-min.js for offline use + optional py3Dmol embed |
| NGLview | MIT | optional `notebook` extra | low | allowed by work order; Jupyter-first users only |
| MolecularNodes | MIT | docs only | low | publication figures in Blender; external tool |
| Avogadro | BSD-3 | docs only | low | manual pose inspection; SDF poses already open natively |

---

## Track 1: flexible receptor & explicit-water docking

### WatVina (`biocheming/watvina`) - RECOMMENDED, ADR-0010

**What it has that DockFlow lacks.**  DockFlow's "hydration" is
hydrogen *addition* (preparator engines); no stage of the pipeline
models displaceable water.  WatVina treats water as a first-class
docking citizen: an FHFT fluid-site occupancy field plus an explicit
rigid-water GCMC network predicts hydration sites from the *dry*
receptor (`--predict_water`), a displacement bias lets those waters
participate in scoring (`--water_bias_weight`), and it outputs an
OpenDX occupancy map plus ranked hydration-site PDBs.  It also brings
pharmacophore generation and template-biased docking (`--genph4` /
`--template`), post-docking MM relaxation, SDF/MOL ligand input (no
PDBQT ligand prep at all), a `--seed` for deterministic reruns, and a
pybind11 API mirroring the CLI.

**What we port.**  Nothing as code - we integrate the *binary* as a
`docking.backend = "watvina"` option in the same pattern as the GNINA
backend (ADR-0001): detect the executable, translate the grid box +
seed + exhaustiveness into CLI flags, parse scores, and record
`backend` + engine version in the manifest.  Two DockFlow-side
additions make it fit: (1) the receptor path passes our prepared
PDBQT *or* the raw PDB (WatVina reads PDB receptors directly - we
record which, because charge models then differ and scores are not
comparable with Vina runs - extend `engine_comparability`), and (2)
water predictions land in `<run>/analysis/hydration_sites.pdb` and the
occupancy map in `<run>/docking/`, surfaced in the report.

**License.**  "Any restrictions to use by non-academics: no license
needed" - free for any use.  Before shipping the download in
`scripts/install_tools.sh`, confirm binary redistribution is welcome
with the author and record the answer in ADR-0010 (the repo stages
prebuilt binaries for Linux/Windows in `bin/`, which suggests yes).

**Why not RxDock first.**  WatVina is self-contained (Eigen, JSON,
force fields compiled in), ships Windows/Linux binaries, and covers
water + pharmacophore + flex side chains in one tool.  RxDock needs a
C++ toolchain build, and its Python tooling is legacy.  If RxDock is
ever added, it is another CLI backend, never an import (LGPL).

### RxDock (`r-dock/rdock`) - DEFER

Empirical scoring with flexible side chains and solvent mapping,
LGPL-2.1, actively maintained by the r-dock org.  DockFlow-relevant
features (cavity-centric pharmacophore restraints, "SF3" scoring) are
attractive, but the value overlaps WatVina's at a higher integration
cost (build from source per platform; no pip; python2-era tooling
around it).  **Port nothing now.**  Revisit if WatVina disappoints in
the water benchmark, or if cavity-definition semantics (rdock's
cavity mapping) are wanted for focused docking.

## Track 2: pharmacophore-guided docking

### PharmacoNet (`SeonghwanSeo/PharmacoNet`) - IDEAS ONLY

Deep-learning protein-based pharmacophore modelling (image instance
segmentation over the pocket) with coarse graph-matching screening -
a *pre-docking filter* that screens millions per CPU-hour, plus a
parameterised pharmacophore-aware scoring function.  MIT, published in
Chem. Sci.  The blocker is the dependency set: PyTorch, Biopython,
`pymol-open-source` (conda-only) and `molvoxel`; a heavy, conda-flavoured
stack that would dominate DockFlow's install size and CI matrix.

**What we take.**  Concepts, not code: (1) the pocket-pointcloud
feature-extraction idea is the right data model for a future
`dockflow pharmacophore` stage; (2) using a fast pharmacophore screen
*before* docking as a library pre-filter fits DockFlow's
virtual-screening workflow (a `--prefilter` idea for the roadmap);
(3) their analytical scoring formula is a candidate consensus
component.  Each of these needs its own ADR if ever built.  Their
sibling GUI (OpenPharmaco) is irrelevant to us.

### OpenPharmacophore (`uibcdf/openpharmacophore`) - VENDOR IF TRACK STARTS

MIT, clean small codebase for pharmacophore *feature perception*
(HBD/HBA/hydrophobic/pos/neg/aromatic) and hypothesis work, built on
RDKit + NGL.  Development stopped mid-2023 (last push 2023-07).  If the
pharmacophore track starts, vendor the feature-perception module (MIT,
with NOTICE attribution) into `dockflow_core/pharmacophore.py` rather
than taking the dependency (dormant upstream + GUI deps we do not
want).  Until then: reference implementation for criteria design.

## Track 3: advanced scoring & plugin architecture

### OpenDock (`guyuehuo/opendock`) - PORT THE PATTERNS

PyTorch docking framework (Bioinformatics 2024, MIT) whose *architecture*
is the valuable part for DockFlow, more than its science: a scoring
**registry** where any `ScoringFunction` can be registered and swapped
into search or rescoring roles, distance-**constraint objects** that
enable covalent/restricted docking by geometry (map- and restraint-based
sampling), and a clean separation between samplers (simulated annealing,
Monte Carlo) and scorers.

**What we port (as MIT-compatible patterns, attribution in the ADR):**
(1) a `ScoringFunction` registry in `ranking.py` so consensus components
(vina, vinardo, gnina_affinity, future pandadock-GNN) register themselves
with a name and a `rescore(poses) -> scores` signature - today's
consensus code is hardcoded to three components; (2) a
`DistanceConstraint` value type for the existing covalent-docking path
(ADR-0004) so bond-length and pocket constraints become data, not
ad-hoc flags; (3) their distance-constraint parsing tests as a model for
ours.  Not ported: the PyTorch samplers/scorers themselves (torch is a
multi-hundred-MB dependency for marginal gain over Vina/GNINA).

### PandaDock (`pritampanda15/PandaDock`) - OPTIONAL ENGINE (LATER)

MIT, pip-installable pure-Python docking suite: MC + quasi-Newton search
with every rotatable bond as an explicit DOF, plus an SE(3)-equivariant
GNN rescoring model trained on 741k complexes.  Two things make it
interesting despite being slower than Vina: it needs **no compiled
engine at all** (a fallback backend where `vina` wheels are missing -
notably Windows/macOS pip installs), and its GNN score is an independent
*opinion* to fold into `ranking.py` consensus.

**What we port.**  Nothing now.  Proposed second step after ADR-0009:
an optional `pandadock` extra + backend adapter (subprocess
`pandadock` CLI or in-process API), plus a `rescore` component for
consensus.  Effort: ~1 day for the adapter + comparability manifest
fields; the scientific validation (does the GNN rank our benchmark
poses sensibly?) is the real cost and should ride on the existing
24-complex redocking benchmark.

## Track 4: 2D ligand-interaction diagrams (new roadmap track)

This track is where the strongest, lowest-risk porting wins live,
because DockFlow already computes interactions (analyzer IFP +
geometric contacts) and only *lacks the rendering layer*.

### ProLIF (`chemosim-lab/ProLIF`) - RECOMMENDED, ADR-0009

**What it has.**  Interaction fingerprints for ligand/protein/DNA/RNA
complexes from docking poses or experimental structures, with a richer
and better-tested interaction vocabulary than DockFlow's homegrown
`ifp_vocabulary` (donor/acceptor-direction hydrogen bonds, hydrophobic,
vdW, salt bridges, cationic/anionic, pi-stacking...), built on RDKit
(no MDAnalysis needed for static structures) - Apache-2.0, active, on
conda-forge and PyPI.  Critically it also ships **LigNetwork**: 2D
interaction diagrams rendered in the browser (D3) - exactly the
"Track 2 roadmap" feature.

**What we port.**  An *adapter, not a rewrite*: when `prolif` is
importable, `analyzer` computes the canonical ProLIF fingerprint for
each top pose (ligand residues vs receptor residues, bit vector per
interacting residue pair) and the report gains a per-pose 2D
LigNetwork diagram (`<run>/visualization/lignetwork_<pose>.html`,
same self-contained-HTML philosophy as the 3Dmol viewer).  The
homegrown IFP stays the always-available fallback; the manifest
records which fingerprint produced `consensus_rank` (they differ in
vocabulary, so ranks are not comparable across engines - the
`comparable_to` machinery gets a fingerprint-flavour entry).  The
consensus IFP term (ranking.py) switches to the ProLIF fingerprint
when present.

**License.**  Apache-2.0 - compatible as an optional dependency;
NOTICE-style attribution in THIRD_PARTY_LICENSES.md.

### PLIP (`pharmai/plip`) - EXTERNAL CLI ONLY, NEVER IMPORT

Rich interaction detection (hydrogen bonds, hydrophobic contacts, salt
bridges, pi-stacking, pi-cation, halogen bonds, water bridges, **metal
coordination**) with PyMOL scene generation and XML/JSON reports - the
reference tool for interaction *reports*.  GPL-2.0: importing plip into
DockFlow would relicence the whole package; **not negotiable**.

**What we port.**  Nothing as code.  Integration route (medium
priority, after Track 4 MVP): run the official Docker image
(`pharmai/plip`) as an optional post-analysis step for users who want
publication-style interaction reports, same process boundary as PyMOL.
The dep-free takeaway we *can* copy legally: their water-bridge and
metal-coordination geometric criteria (already listed in
`docs/interaction_criteria.md` conventions table) sharpen our own
heuristics - criteria/facts are not copyrightable, code is.

## Track 5: 3D visualisation & rendering

### 3Dmol.js / py3Dmol (`3dmol/3Dmol.js`) - ALREADY USED; TWO UPGRADES

DockFlow's `interactive_viewer.py` already renders with 3Dmol.js from
CDN (ADR-0006).  Two cheap improvements fall out of this review, both
BSD-3 and inside the work order's allowed dependency set:
(1) **vendor `3Dmol-min.js`** into `docs/assets/` and emit a
file-relative `<script src>` with the CDN version as fallback - today's
page breaks behind firewalls, contradicting the "open offline" claim in
its own footer; (2) **py3Dmol embed for notebook users** - a tiny
`notebook_view(run_dir)` helper (optional `notebook` extra) that shows
receptor + poses inline in Jupyter.  Effort: hours, not days.

### NGLview (`nglviewer/nglview`) - OPTIONAL, JOINTLY WITH py3Dmol

MIT Jupyter widget (NGL renderer) with first-class MDAnalysis/RDKit
integration.  The work order pre-approved NGL as a dependency.
Decision: offer **either** py3Dmol (lighter, embeds in HTML) **or**
nglview (richer, needs a running kernel) in the same optional
`notebook` extra; default to py3Dmol so the no-Jupyter path never
regresses.  No code ported.

### Molecular Nodes (`BradyAJohnston/MolecularNodes`) - DOCS ONLY

MIT Blender add-on for publication-grade molecular rendering and
animation (now powered by Biotite + MDAnalysis).  DockFlow's SDF pose
export (ADR-0002) already feeds it: `import_pdb`/MDAnalysis reads our
receptor + pose SDFs directly.  **Action:** a short "Publication
figures with Molecular Nodes" cookbook section in USER_GUIDE.md
pointing at the exported files.  No dependency, no porting.

### Avogadro 2 (`OpenChemistry/avogadro`) - DOCS ONLY

BSD-3 C++/Qt molecular editor - the right answer for *manual* pose
inspection/editing and quick PDBQT/SDF eyeballing without PyMOL.
DockFlow's outputs (PDBQT, exported pose SDFs) open natively.
**Action:** mention in USER_GUIDE troubleshooting as the lightweight
desktop viewer option (PyMOL install avoided).  No dependency, no
porting.

---

## Recommended sequence (and what it costs)

1. **Now (0.3.x):** vendor `3Dmol-min.js` for offline interactive.html
   + offline footer fix.  Hours.  No ADR needed (allowed dependency).
2. **0.4 - Track 2 MVP:** ADR-0009 ProLIF optional backend + LigNetwork
   per-pose diagrams + fingerprint provenance in manifest.  ~2-3 days.
3. **0.4 - water track spike:** ADR-0010 WatVina backend behind
   `--backend watvina`; run the 24-complex benchmark with water
   predictions on; publish comparability caveats.  ~2 days + validation
   compute.
4. **0.5:** OpenDock pattern port (scoring registry + constraint
   objects) while refactoring ranking.py; PandaDock optional engine +
   GNN consensus component (both gated on ADRs + benchmark evidence).
5. **Deferred:** PLIP docker report step; RxDock backend; pharmacophore
   stage (ideas from PharmacoNet; vendor OpenPharmacophore perception
   module if it starts).

Each step keeps the standing principle intact: mechanical decisions
automated, scientific decisions (which engine, which fingerprint,
which waters to trust) exposed in the manifest, recorded, and
overridable.
