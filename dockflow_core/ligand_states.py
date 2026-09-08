"""Ligand state enumeration: protonation, tautomers, stereochemistry, salts.

Work order items 11-14.  Every function here is *state chemistry with
explicit provenance*: it either returns the states it produced (to be
docked individually and recorded in the manifest) or a report about the
input state, so no chemical assumption is ever silent.

* protonation (item 11): dimorphite-dl at a configured pH; ``major``
  mode records the dominant microspecies, ``enumerate`` mode docks every
  state dimorphite-dl considers populated at that pH.
* tautomers (item 12): RDKit's ``TautomerEnumerator`` (canonical set;
  no pH weighting is applied by RDKit - the limitation is recorded).
* stereocenters (item 13): undefined chiral centers are detected and
  reported; optional enumeration docks every stereoisomer.
* salts (item 14): fragment analysis records exactly what was stripped.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

_HALIDE_LIKE_MAX_HEAVY = 2  # F/Cl/Br/I counter-ions and simple anions


@dataclass
class StateVariant:
    """One dockable state of a ligand (protonation/tautomer/stereo level)."""

    mol: Any
    tag: str = ""            # suffix for the ligand id, e.g. "prot1"
    kind: str = "input"      # input | protonation | tautomer | stereoisomer
    index: int = 1
    smiles: str = ""
    detail: dict[str, Any] = field(default_factory=dict)

    def describe(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "index": self.index,
            "tag": self.tag,
            "smiles": self.smiles,
            **self.detail,
        }


# ---------------------------------------------------------------------------
# Item 11: protonation (dimorphite-dl)
# ---------------------------------------------------------------------------
def protonation_states(
    mol: Any,
    mode: str = "major",
    ph: float = 7.4,
    window: float = 0.5,
    max_variants: int = 8,
) -> tuple[list[Any], dict[str, Any]]:
    """Protonation states of ``mol`` at pH ``ph`` +/- ``window``.

    ``mode``: ``"major"`` (record the dominant microspecies, dock one
    state) or ``"enumerate"`` (dock every populated state).

    Returns ``(mols, report)`` where ``mols`` is the list of RDKit
    molecules to dock (the input molecule unchanged when dimorphite-dl
    is unavailable) and ``report`` is the provenance dict for the
    manifest.  Failures are recorded, never raised: the input state is
    used and the manifest says so.
    """
    from rdkit import Chem

    report: dict[str, Any] = {
        "mode": mode,
        "ph": ph,
        "window": window,
        "source": "input protonation state (no pH adjustment)",
        "major_microspecies_smiles": None,
        "net_charge": None,
        "states": [],
    }
    input_smiles = Chem.MolToSmiles(mol)
    report["states"].append({"smiles": input_smiles, "source": "input",
                             "index": 1})
    try:
        from dimorphite_dl import protonate_smiles
    except ImportError:
        report["source"] = ("input protonation state "
                            "(dimorphite-dl not installed)")
        report["note"] = "install dimorphite-dl for pH-aware protonation"
        return [mol], report

    try:
        variants = protonate_smiles(
            input_smiles, ph_min=ph - window, ph_max=ph + window,
            max_variants=max_variants,
        )
        variants = [v for v in variants if v]
    except Exception as exc:  # noqa: BLE001 - optional chemistry, never fatal
        report["source"] = f"input state (protonation failed: {exc})"
        return [mol], report

    if not variants:
        report["source"] = "input state (dimorphite-dl returned nothing usable)"
        return [mol], report

    mols = []
    states = []
    for index, smiles in enumerate(variants[:max_variants], start=1):
        variant_mol = Chem.MolFromSmiles(smiles)
        if variant_mol is None:
            continue
        mols.append(variant_mol)
        states.append({
            "smiles": smiles,
            "index": index,
            "net_charge": Chem.GetFormalCharge(variant_mol),
            "source": "dimorphite-dl",
        })
    if not mols:
        report["source"] = "input state (no usable dimorphite-dl output)"
        return [mol], report

    report["source"] = (f"dimorphite-dl pH {ph} +/- {window} "
                        f"({'major microspecies' if mode == 'major' else 'enumeration'})")
    report["states"] = states
    report["major_microspecies_smiles"] = states[0]["smiles"]
    report["net_charge"] = states[0]["net_charge"]
    if mode == "major":
        if len(mols) > 1:
            report["note"] = (
                f"{len(mols)} states populated at this pH; docking the "
                f"major microspecies only (use protonate: enumerate to "
                f"dock all)"
            )
        return [mols[0]], report
    return mols, report


# ---------------------------------------------------------------------------
# Item 12: tautomers (RDKit TautomerEnumerator)
# ---------------------------------------------------------------------------
def tautomer_states(
    mol: Any,
    max_tautomers: int = 16,
) -> tuple[list[Any], dict[str, Any]]:
    """Enumerate the canonical tautomer set of ``mol``.

    RDKit's ``TautomerEnumerator`` produces the structurally possible
    tautomers without pH weighting; the returned report records that
    limitation.  The input tautomer is always included (index 1).
    """
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    report: dict[str, Any] = {
        "source": "single input tautomer (no enumeration)",
        "tautomers": [],
        "note": None,
    }
    try:
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(max_tautomers)
        tautomers = list(enumerator.Enumerate(Chem.Mol(mol)))
    except Exception as exc:  # noqa: BLE001
        report["note"] = f"tautomer enumeration failed: {exc}"
        return [mol], report

    tautomers = tautomers[:max_tautomers]
    if len(tautomers) <= 1:
        report["tautomers"] = [
            {"index": 1, "smiles": Chem.MolToSmiles(mol), "source": "input"}]
        return [mol], report

    canonical_input = Chem.MolToSmiles(mol)
    ordered = [t for t in tautomers if Chem.MolToSmiles(t) != canonical_input]
    input_first = [Chem.Mol(mol)] + ordered
    report["source"] = (f"RDKit TautomerEnumerator, canonical set "
                        f"({len(input_first)} tautomers, capped at "
                        f"{max_tautomers}; no pH weighting - RDKit "
                        f"enumerates structural possibilities)")
    report["tautomers"] = [
        {"index": i, "smiles": Chem.MolToSmiles(t),
         "source": "input" if i == 1 else "enumerated"}
        for i, t in enumerate(input_first, start=1)
    ]
    return input_first, report


# ---------------------------------------------------------------------------
# Item 13: stereocenters
# ---------------------------------------------------------------------------
def stereocenter_report(mol: Any) -> dict[str, Any]:
    """Defined/undefined chiral centers + stereoisomer count (item 13a)."""
    from rdkit import Chem

    try:
        centers = Chem.FindMolChiralCenters(
            mol, useLegacyImplementation=False, includeUnassigned=True)
    except Exception:  # noqa: BLE001 - very old RDKit
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    undefined = [(index, tag) for index, tag in centers if tag == "?"]
    possible = 2 ** len(undefined) if undefined else 0
    return {
        "defined": len(centers) - len(undefined),
        "undefined": len(undefined),
        "undefined_centers": [index for index, _ in undefined],
        "possible_stereoisomers": possible,
    }


def stereoisomer_states(
    mol: Any,
    max_isomers: int = 8,
) -> tuple[list[Any], dict[str, Any]]:
    """Enumerate stereoisomers for the unassigned centers (item 13c)."""
    from rdkit import Chem
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

    report = stereocenter_report(mol)
    if not report["undefined"]:
        report["source"] = "no undefined stereocenters"
        return [mol], report
    options = StereoEnumerationOptions(
        onlyUnassigned=True, unique=True, maxIsomers=max_isomers)
    isomers = list(EnumerateStereoisomers(mol, options=options)) or [mol]
    isomers = isomers[:max_isomers]
    report["source"] = (f"RDKit EnumerateStereoisomers (onlyUnassigned, "
                        f"unique; {len(isomers)} of "
                        f"{report['possible_stereoisomers']} possible, "
                        f"capped at {max_isomers})")
    report["enumerated"] = [
        {"index": i, "smiles": Chem.MolToSmiles(iso)}
        for i, iso in enumerate(isomers, start=1)
    ]
    return isomers, report


# ---------------------------------------------------------------------------
# Item 14: salt / fragment analysis
# ---------------------------------------------------------------------------
def salt_report(mol: Any) -> dict[str, Any] | None:
    """Fragment analysis for salt stripping provenance (item 14a).

    Returns ``None`` for single-fragment molecules.  The kept fragment is
    the largest (matching the preparation behaviour); removed fragments
    are listed with SMILES and heavy-atom counts.  A
    ``large_fragment_warning`` is included when a removed fragment is
    bigger than a typical counter-ion (halide / small anion).
    """
    from rdkit import Chem

    fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
    if len(fragments) <= 1:
        return None
    ranked = sorted(fragments, key=lambda m: m.GetNumAtoms(), reverse=True)
    kept = ranked[0]
    removed = ranked[1:]
    removed_info = []
    for fragment in removed:
        try:
            smiles = Chem.MolToSmiles(fragment)
        except Exception:  # noqa: BLE001 - fragment may fail sanitisation
            smiles = f"<{fragment.GetNumAtoms()} atoms>"
        removed_info.append({
            "smiles": smiles,
            "heavy_atoms": fragment.GetNumHeavyAtoms(),
        })
    large = [info for info in removed_info
             if info["heavy_atoms"] > _HALIDE_LIKE_MAX_HEAVY]
    report = {
        "removed_fragments": removed_info,
        "kept_fragment_atoms": kept.GetNumAtoms(),
        "kept_fragment_heavy_atoms": kept.GetNumHeavyAtoms(),
        "kept_fragment_smiles": Chem.MolToSmiles(kept),
    }
    if large:
        names = ", ".join(
            f"{info['smiles']} ({info['heavy_atoms']} heavy atoms)"
            for info in large[:3]
        )
        report["large_fragment_warning"] = (
            f"Salt fragment(s) {names} stripped; confirm they are "
            f"counter-ions, not a second pharmacophore."
        )
    return report
