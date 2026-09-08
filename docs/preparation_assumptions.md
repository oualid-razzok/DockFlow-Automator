# Preparation assumptions

Every preparation decision DockFlow makes on the way from a deposited
PDB entry (or SDF/SMILES ligand) to a dockable PDBQT.  Each entry
follows the same three-sentence structure the final review asked for:
**what DockFlow does**, **what it does NOT do**, and **what the
researcher should do for publication-grade work**.

The machine-readable record of the choices actually made in a run is
`manifest.json` → `receptor.decisions` and `ligands[*].prep`; the
report footer ("Assumptions this run made") summarises them, and a
catalytic-residue warning fires when pH-sensitive residues sit near the
binding site (see `dockflow_core/structure_qc.py`).

The guiding principle: DockFlow automates the *mechanical* decisions and
records them, but the *scientific* decisions (protonation, tautomer,
stereochemistry, keeping cofactors) remain the researcher's to override.

## Receptor preparation

### Hydrogen addition

* **Does:** keeps the receptor's DEPOSITED protonation states as solved
  and adds hydrogens by toolkit template at a nominal pH 7.4 (OpenBabel
  or RDKit residue-template lookup; the engine actually used is recorded
  in `manifest.receptor.decisions`).
* **Does NOT:** run a pKa calculation.  Crystallographic protonation is
  rarely resolved at typical PDB resolution; "pH 7.4 nominal" describes
  the template table, not the local electrostatic environment.
  Residues whose pKa is shifted by the active site (catalytic Asp/Glu,
  His tautomers, buried Lys) may therefore be mis-protonated.
* **Should:** run PROPKA or H++ on the receptor, inspect the
  active-site residues, and pass the pKa-informed (re)protonated PDB as
  the DockFlow input.  The report's catalytic-residue warning names the
  His/Asp/Glu residues near the box that deserve that check.

### Tautomers (receptor)

* **Does:** nothing — protein receptors are docked as the single
  deposited/templated protonation state; His is typed by whatever
  tautomer the hydrogen addition produced.
* **Does NOT:** enumerate or score His tautomers.
* **Should:** for catalytic His residues, test both tautomers manually
  (two runs, one per tautomer) when the tautomer matters.

### Tautomers (ligand)

* **Does:** docks the SINGLE input tautomer (Meeko receives exactly what
  the input file/SMILES contains).
* **Does NOT:** enumerate tautomers unless explicitly configured
  (ligand preparation records "single input tautomer (no enumeration)"
  in `manifest.ligands[*].prep.tautomer`).
* **Should:** for ligands with tautomeric groups (imidazole, tetrazole,
  amide/phenol pairs), dock each plausible tautomer as a separate ligand
  and compare, or use a dedicated tautomer-enumeration tool upstream.

### Stereochemistry

* **Does:** docks ligands exactly as input; defined stereocenters are
  preserved (recorded as `stereocenters.defined` /
  `possible_stereoisomers` in the ligand prep record).
* **Does NOT:** enumerate undefined stereocenters unless configured
  (`enumerate_stereoisomers`); undefined centers raise a WARNING, never
  a silent guess.
* **Should:** provide fully specified stereochemistry (isomeric SMILES
  or a 3D SDF) — racemic/undefined input means the docked ensemble
  represents an unspecified mixture.

### Charge model

* **Does:** assigns Gasteiger partial charges (OpenBabel engine) or
  Gasteiger + formal charges (RDKit engine); the model is recorded in
  `manifest.receptor.decisions.charge_model`.
* **Does NOT:** derive charges from quantum mechanics or from
  RESP/AM1-BESP fitting; metal sites, heme groups and unusual oxidation
  states get placeholder-quality charges.
* **Should:** for metal-dependent or highly polarised sites, rescore
  top poses with a QM-derived charge model or a scoring function with
  metal terms (AD4 metal maps, GNINA-CNN).

### Waters

* **Does:** removes all water molecules by default and records the count
  (`manifest.receptor.decisions.waters`).
* **Does NOT:** decide which waters are conserved/structural — a
  bridging water can be essential for the pose (it is removed with
  everything else).
* **Should:** inspect conserved-water analysis (e.g. WaterMap-style
  tools, or literature for the target) and, if a bridge water matters,
  re-add it explicitly to the receptor input or use a water-aware
  docking backend (see ROADMAP_POST_V1.md for the WatVina plan).

### Metals and cofactors

* **Does:** removes HETATM species (metals, cofactors, ions) by default,
  recording each removed species and count in
  `manifest.receptor.decisions.metals_cofactors`; kept metals get a
  coordination-environment record (`metal_coordination`) and a report
  warning.
* **Does NOT:** know which cofactor is catalytically required.  A
  structural Zn²⁺ removed from a metalloprotease changes the pocket.
* **Should:** pass known functional metals/cofactors via
  `receptor.keep_resnames` (e.g. `ZN,HEME,MG`) so they are kept and
  warned about, and prefer metal-aware scoring for such sites.

### Alternate locations (altLoc)

* **Does:** resolves alternate conformations to the highest-occupancy
  variant by default (`altloc: best`); other variants are dropped.
* **Does NOT:** model partial occupancy or ensemble docking of the
  alternate conformers.
* **Should:** pass `altloc: "A"` (or another specific altLoc) when the
  deposited alternate is the biologically relevant one; the policy used
  is recorded in `manifest.receptor.decisions.alternate_conformations`.

### Disulfides, missing residues, chain gaps

* **Does:** detects disulfides (S–S < 2.5 Å), flags reduced cysteines
  and numbering gaps (missing residues) heuristically, and records them
  in the manifest with report warnings.
* **Does NOT:** repair anything — Vina cannot model loops.
* **Should:** rebuild long gaps with Modeller/SWISS-MODEL before
  docking if the gap borders the binding site (the warning says so).

## Ligand preparation (Meeko)

* **Does:** converts the input molecule to a PDBQT with Meeko (AD4 atom
  types, Gasteiger charges, torsion tree with macrocycle breaking);
  salts are stripped by Meeko default; embedding/minimisation (when
  input is SMILES) is seeded and recorded.
* **Does NOT:** adjust the protonation state (input state is kept; the
  optional dimorphite-dl integration must be configured explicitly and
  is a *choice*, not a default), enumerate tautomers or stereoisomers
  (see above), or keep co-salts.
* **Should:** curate the input protonation/tautomer/stereo state before
  docking — for publication-grade work this is the single highest-impact
  ligand-side decision, and it is deliberately left to the researcher.

## Grid box

* **Does:** derives the search box from the co-crystallized ligand of
  the target (4 Å padding) when one exists; the box's provenance
  (source, padding, assumption strength) is recorded in
  `manifest.gridbox` and the report.
* **Does NOT:** know the pocket when no reference ligand exists — the
  whole-structure fallback is labelled with the weakest assumption
  strength and warned about.
* **Should:** provide explicit `gridbox.center/size` or a curated
  active-site residue list when the pocket is known from the literature;
  confirm weak/medium boxes with conservation analysis.

## Where the assumptions are surfaced

| Surface | Content |
|---|---|
| `report.md` | assumptions footer, catalytic-residue warning, grid-box assumption strength, degraded-mode banner |
| `manifest.json` | full machine-readable decisions per stage |
| `docs/interaction_criteria.md` | thresholds and their convention/heuristic classification |
| `benchmarks/` | the measured consequences of these assumptions (redocking, MGLTools comparison, enrichment) |
