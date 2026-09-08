"""Receptor structure quality checks (work order items 6, 7 and 10).

These checks are *heuristic structure checks*, not simulations:

* missing-residue detection - numbering gaps in the ATOM records (no
  SEQRES cross-check; author numbering can legally skip numbers),
* disulfide detection - SG-SG distance < 2.5 A between cysteines,
* reduced-cysteine flagging - free Cys SG near another free Cys SG
  (proximity heuristic, no redox context),
* metal-coordination detection - kept metal centres with their
  coordinating residues within 2.8 A and a coordination-geometry label.

Every finding is recorded in the run manifest
(``manifest.receptor.decisions``) and surfaced as a report warning so
gapped / reduced / metal-dependent receptors are never docked silently.
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Sequence

    from .pdbio import Atom

# Standard amino acids (three-letter codes) - used to walk polymer gaps.
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",  # selenomethionine reads as a polymer residue
    "SEC", "PYL",
}

WATER_RESNAMES = {"HOH", "DOD", "WAT", "H2O", "TIP", "TIP3", "SOL", "W"}

# Metal elements commonly found in protein active sites.  Matched against
# the element column first, the residue name second (covers names like
# FE2 / CU1 /ZN missing element columns).
METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "TI", "V", "CR",
    "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG", "AL", "GA", "IN", "SN",
    "PB", "SB", "BI", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG",
    "LA", "CE", "ND", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
    "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "TH", "U",
}
METAL_RESNAMES = set(METAL_ELEMENTS) | {
    "FE2", "FE3", "CU1", "CU2", "MN3", "ZN2", "CO2", "NI2", "CD2", "HG2",
}


def _is_metal(atom: Atom) -> bool:
    element = (atom.element or "").strip().upper()
    if element:
        return element in METAL_ELEMENTS
    return (atom.resname or "").strip().upper() in METAL_RESNAMES


def _is_polymer(atom: Atom) -> bool:
    return ((atom.record_type or "ATOM").strip().upper() == "ATOM"
            and (atom.resname or "").strip().upper() in STANDARD_RESIDUES)


def _xyz(atoms: Sequence[Atom]) -> np.ndarray:
    return np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)


# ---------------------------------------------------------------------------
# Item 6: missing residues
# ---------------------------------------------------------------------------
def detect_missing_residues(atoms: Sequence[Atom],
                            max_gap_to_report: int = 30) -> list[dict]:
    """Numbering gaps in the polymer of each chain.

    A *gap* is two consecutive polymer residues (sorted by resseq, then
    insertion code) whose numbering jumps, with nothing in between.  This
    is a heuristic: crystallographers may skip numbers without unmodelled
    residues, and terminal disorder is invisible to it.  Gaps larger than
    ``max_gap_to_report`` are reported as one entry with an approximate
    count (huge numbering jumps usually mean separate domains or author
    numbering conventions, not 200 missing residues).

    Returns ``[{chain, resseq_range, count}, ...]``.
    """
    per_chain: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for atom in atoms:
        if _is_polymer(atom):
            per_chain[atom.chain.strip() or "_"].append(
                (int(atom.resseq), (atom.icode or "").strip())
            )
    gaps: list[dict] = []
    for chain, residues in per_chain.items():
        keys = sorted(set(residues))
        for (resseq, icode), (next_seq, next_icode) in zip(keys, keys[1:], strict=False):
            same_seq_next_icode = resseq == next_seq and icode < next_icode
            contiguous = (resseq == next_seq - 1 and next_icode == "") or \
                (resseq == next_seq - 1 and next_icode != "" and icode == "") or \
                same_seq_next_icode
            if contiguous:
                continue
            missing = next_seq - resseq - 1
            if missing <= 0:
                continue
            first, last = resseq + 1, next_seq - 1
            span = str(first) if first == last else f"{first}-{last}"
            if missing > max_gap_to_report:
                gaps.append({
                    "chain": chain,
                    "resseq_range": span,
                    "count": missing,
                    "note": "large numbering jump (approximate; may be "
                            "author numbering, not unmodelled residues)",
                })
            else:
                gaps.append({
                    "chain": chain,
                    "resseq_range": span,
                    "count": missing,
                })
    return gaps


def missing_residue_warning(gaps: Sequence[dict]) -> str | None:
    """Human-readable warning for report.md (peer item 6b)."""
    if not gaps:
        return None
    total = sum(int(gap["count"]) for gap in gaps)
    detail = "; ".join(
        f"chain {gap['chain']}: {gap['resseq_range']}" for gap in gaps[:4]
    )
    more = f" (+{len(gaps) - 4} more)" if len(gaps) > 4 else ""
    return (
        f"Receptor has {total} missing residues ({detail}{more}). "
        "Vina cannot model loops; consider Modeller/SWISS-MODEL "
        "before docking. (Heuristic numbering-gap detection; see "
        "manifest.receptor.decisions.missing_residues.)"
    )


# ---------------------------------------------------------------------------
# Item 7: disulfides
# ---------------------------------------------------------------------------
def detect_disulfides(atoms: Sequence[Atom],
                      cutoff: float = 2.5) -> list[dict]:
    """Detect disulfide bridges: CYS SG-SG pairs closer than ``cutoff`` A.

    Returns ``[{chain_a, resseq_a, chain_b, resseq_b, distance_a}, ...]``
    with each bridge listed once (a-b sorted).
    """
    sulfurs = [
        atom for atom in atoms
        if (atom.resname or "").strip().upper() == "CYS"
        and (atom.name or "").strip().upper() == "SG"
    ]
    bridges: list[dict] = []
    for i, first in enumerate(sulfurs):
        for second in sulfurs[i + 1:]:
            if first is second:
                continue
            same_residue = (first.chain == second.chain
                            and first.resseq == second.resseq)
            if same_residue:
                continue
            distance = float(np.linalg.norm(
                np.array([first.x, first.y, first.z])
                - np.array([second.x, second.y, second.z])))
            if distance <= cutoff:
                a, b = sorted(
                    [(first.chain.strip(), int(first.resseq)),
                     (second.chain.strip(), int(second.resseq))])
                bridges.append({
                    "chain_a": a[0], "resseq_a": a[1],
                    "chain_b": b[0], "resseq_b": b[1],
                    "distance_a": round(distance, 2),
                })
    return bridges


def disulfide_partner_map(disulfides: Sequence[dict]) -> dict[tuple, tuple]:
    """Map (chain, resseq) -> partner (chain, resseq) over all bridges."""
    partners: dict[tuple, tuple] = {}
    for bridge in disulfides:
        a = (bridge["chain_a"], int(bridge["resseq_a"]))
        b = (bridge["chain_b"], int(bridge["resseq_b"]))
        partners[a] = b
        partners[b] = a
    return partners


def detect_reduced_cysteines(atoms: Sequence[Atom],
                             disulfides: Sequence[dict],
                             proximity: float = 8.0) -> list[dict]:
    """Free cysteines whose SG sits near another free cysteine SG.

    Heuristic (peer item 7c): a free Cys near another free Cys in an
    oxidising context *may* be a missed disulfide or need protonation
    state review.  No redox environment is modelled - the finding is a
    flag for manual review, not a correction.
    """
    partners = disulfide_partner_map(disulfides)
    sulfurs = [
        atom for atom in atoms
        if (atom.resname or "").strip().upper() == "CYS"
        and (atom.name or "").strip().upper() == "SG"
        and (atom.chain.strip(), int(atom.resseq)) not in partners
    ]
    flagged: list[dict] = []
    seen: set[tuple] = set()
    for i, first in enumerate(sulfurs):
        for second in sulfurs[i + 1:]:
            key = tuple(sorted([
                (first.chain.strip(), int(first.resseq)),
                (second.chain.strip(), int(second.resseq))]))
            if key in seen:
                continue
            distance = float(np.linalg.norm(
                np.array([first.x, first.y, first.z])
                - np.array([second.x, second.y, second.z])))
            if distance <= proximity:
                seen.add(key)
                flagged.append({
                    "chain": first.chain.strip(),
                    "resseq": int(first.resseq),
                    "near_chain": second.chain.strip(),
                    "near_resseq": int(second.resseq),
                    "sg_sg_distance_a": round(distance, 2),
                })
    return flagged


def reduced_cys_warning(flagged: Sequence[dict]) -> str | None:
    """Human-readable warning for report.md (peer item 7c)."""
    if not flagged:
        return None
    names = ", ".join(
        f"Cys{row['resseq']}" +
        (f"/Cys{row['near_resseq']}" if row["near_resseq"] != row["resseq"]
         else "")
        for row in flagged[:4]
    )
    return (
        f"{names} appear(s) reduced (free SG, no S-S partner, another free "
        f"Cys within 8 A); may need protonation state review. "
        "(Heuristic: see manifest.receptor.decisions.reduced_cysteines.)"
    )


def broken_disulfide_warning(atoms_out: Sequence[Atom],
                             disulfides: Sequence[dict]) -> str | None:
    """Warn when a detected S-S bridge lost an SG atom in preparation."""
    surviving = {
        (atom.chain.strip(), int(atom.resseq))
        for atom in atoms_out
        if (atom.resname or "").strip().upper() == "CYS"
        and (atom.name or "").strip().upper() == "SG"
    }
    broken = [
        bridge for bridge in disulfides
        if (bridge["chain_a"], int(bridge["resseq_a"])) not in surviving
        or (bridge["chain_b"], int(bridge["resseq_b"])) not in surviving
    ]
    if not broken:
        return None
    detail = ", ".join(
        f"{b['chain_a']}{b['resseq_a']}-{b['chain_b']}{b['resseq_b']}"
        for b in broken[:4]
    )
    return (
        f"Disulfide bridge(s) {detail} lost an SG atom during preparation; "
        "the receptor geometry changed relative to the crystal structure."
    )


# ---------------------------------------------------------------------------
# Item 10: metal coordination
# ---------------------------------------------------------------------------
def detect_metal_coordination(atoms: Sequence[Atom],
                              cutoff: float = 2.8) -> list[dict]:
    """Coordination environment of every retained metal atom.

    Coordinating residues are protein/cofactor atoms within ``cutoff`` A
    of the metal centre (waters excluded - they are removed by default
    anyway).  The geometry label is a coordination-number heuristic:
    6 -> octahedral, 4 -> tetrahedral (square planar is not
    distinguished), anything else -> "other".

    Returns ``[{metal_resname, chain, resseq, coordinating_residues,
    coordination_number, geometry_heuristic}, ...]``.
    """
    metals = [atom for atom in atoms if _is_metal(atom)
              and (atom.resname or "").strip().upper() not in WATER_RESNAMES]
    if not metals:
        return []
    others = [
        atom for atom in atoms
        if atom not in metals
        and (atom.resname or "").strip().upper() not in WATER_RESNAMES
    ]
    other_xyz = _xyz(others)
    findings: list[dict] = []
    for metal in metals:
        center = np.array([metal.x, metal.y, metal.z])
        if not len(other_xyz):
            break
        distances = np.linalg.norm(other_xyz - center, axis=1)
        ligand_indices = np.where(distances <= cutoff)[0]
        residues: list[dict] = []
        seen: set[tuple] = set()
        for index in ligand_indices:
            atom = others[int(index)]
            key = (atom.chain.strip(), int(atom.resseq),
                   (atom.resname or "").strip().upper())
            if key not in seen:
                seen.add(key)
                residues.append({
                    "chain": key[0], "resseq": key[1], "resname": key[2],
                })
        number = len(ligand_indices)
        if number == 6:
            geometry = "octahedral"
        elif number == 4:
            geometry = "tetrahedral"
        else:
            geometry = "other"
        findings.append({
            "metal_resname": (metal.resname or "").strip().upper(),
            "chain": metal.chain.strip(),
            "resseq": int(metal.resseq),
            "coordinating_residues": residues,
            "coordination_number": number,
            "geometry_heuristic": geometry,
            "coordinating_within_a": cutoff,
        })
    return findings


def metal_warning(findings: Sequence[dict],
                  backend: str | None) -> str | None:
    """Backend-aware metal warning for report.md (peer item 10b)."""
    if not findings:
        return None
    names = ", ".join(sorted({
        f"{row['metal_resname']} ({row['geometry_heuristic']}, "
        f"{row['coordination_number']} ligands)"
        for row in findings
    }))
    message = (
        f"Metal site(s) retained: {names}. Coordination environment "
        "recorded in manifest.receptor.decisions.metal_coordination "
        "(distance heuristic only)."
    )
    if backend in (None, "", "auto", "vina", "python"):
        message += (
            " Vina has no metal-coordination scoring term; consider "
            "GNINA-CNN or AD4 metal maps for metal-dependent sites."
        )
    return message


# ---------------------------------------------------------------------------
# Catalytic / pH-sensitive residues near the binding pocket
# (final tightening pass, item 6: protonation honesty)
# ---------------------------------------------------------------------------
# Residues whose protonation state is most often wrong after template-based
# hydrogen addition: His (two tautomers, pKa ~6.5), catalytic Asp/Glu
# (pKa shifts of several units in active sites) and the Ser-His-Asp/Glu
# catalytic triad of hydrolases/proteases.
_HIS_RING_N = ("ND1", "NE2")
_ACIDIC_O = ("OD1", "OD2", "OE1", "OE2")
_TRIAD_SER_OG_TO_HIS_N_A = 4.0
_TRIAD_ACID_O_TO_HIS_N_A = 4.5


def _residue_atoms(atoms: Sequence[Atom]) -> dict[tuple[str, int, str], list[Atom]]:
    """Group standard-residue atoms by (chain, resseq, resname)."""
    groups: dict[tuple[str, int, str], list[Atom]] = defaultdict(list)
    for atom in atoms:
        resname = atom.resname.strip().upper()
        if resname in STANDARD_RESIDUES and atom.chain:
            groups[(atom.chain.strip(), int(atom.resseq), resname)].append(atom)
    return groups


def detect_catalytic_residues(
    atoms: Sequence[Atom],
    box_center: Sequence[float] | None = None,
    box_size: Sequence[float] | None = None,
    pocket_shell: float = 6.0,
) -> list[dict]:
    """Heuristic detection of pH-sensitive / catalytic residues.

    Geometry-only heuristics (no catalytic-residue database lookup):

    * **catalytic triad**: a Ser-His-Asp/Glu arrangement anywhere in the
      receptor - Ser OG within 4.0 A of a His ring nitrogen, and an
      Asp/Glu carboxylate oxygen within 4.5 A of a His ring nitrogen
      (the classic hydrolase/protease charge-relay geometry);
    * **pocket pH-sensitive residues**: His or Asp/Glu side chains with
      any atom inside the grid box or within ``pocket_shell`` A of it
      (requires ``box_center``/``box_size``).

    These are exactly the residues whose protonation state is most likely
    wrong after template-based hydrogen addition (pH 7.4 nominal, NOT a
    pKa calculation) - see docs/preparation_assumptions.md.
    """
    groups = _residue_atoms(atoms)
    findings: list[dict] = []

    # ---- catalytic triad (Ser-His-Asp/Glu charge relay) ------------------
    his_residues = [(key, atoms_) for key, atoms_ in groups.items()
                    if key[2] == "HIS"]
    ser_residues = [atom for key, atoms_ in groups.items()
                    if key[2] == "SER"
                    for atom in atoms_ if atom.name.strip() == "OG"]
    acid_oxygens = [atom for key, atoms_ in groups.items()
                    if key[2] in ("ASP", "GLU")
                    for atom in atoms_ if atom.name.strip() in _ACIDIC_O]
    for his_key, his_atoms in his_residues:
        his_n = [a for a in his_atoms if a.name.strip() in _HIS_RING_N]
        if not his_n:
            continue
        near_ser = [
            ser for ser in ser_residues
            if min(np.linalg.norm([ser.x - n.x, ser.y - n.y, ser.z - n.z])
                   for n in his_n) <= _TRIAD_SER_OG_TO_HIS_N_A]
        if not near_ser:
            continue
        near_acid = [
            oxy for oxy in acid_oxygens
            if min(np.linalg.norm([oxy.x - n.x, oxy.y - n.y, oxy.z - n.z])
                   for n in his_n) <= _TRIAD_ACID_O_TO_HIS_N_A]
        if not near_acid:
            continue
        ser_key = (near_ser[0].chain.strip(), int(near_ser[0].resseq),
                   near_ser[0].resname.strip().upper())
        oxy = near_acid[0]
        acid_key = (oxy.chain.strip(), int(oxy.resseq),
                    oxy.resname.strip().upper())
        findings.append({
            "type": "catalytic_triad",
            "residues": [
                f"Ser {ser_key[0]}{ser_key[1]}",
                f"His {his_key[0]}{his_key[1]}",
                f"{acid_key[2]} {acid_key[0]}{acid_key[1]}",
            ],
            "detail": "Ser-His-Asp/Glu charge-relay geometry detected",
        })

    # ---- pH-sensitive residues near the pocket ---------------------------
    if box_center is not None and box_size is not None:
        center = np.asarray(box_center, dtype=float)
        half = np.asarray(box_size, dtype=float) / 2.0 + pocket_shell
        lower, upper = center - half, center + half
        for key, atoms_ in groups.items():
            if key[2] not in ("HIS", "ASP", "GLU"):
                continue
            inside = [
                a for a in atoms_
                if (lower[0] <= a.x <= upper[0]
                    and lower[1] <= a.y <= upper[1]
                    and lower[2] <= a.z <= upper[2])
            ]
            if inside:
                findings.append({
                    "type": "pocket_ph_sensitive",
                    "residues": [f"{key[2]} {key[0]}{key[1]}"],
                    "detail": (f"{key[2]} side chain within "
                               f"{pocket_shell:.0f} A of the grid box"),
                })
    return findings


def catalytic_residue_warning(findings: Sequence[dict]) -> str | None:
    """Protonation-honesty warning for report.md (final pass, item 6)."""
    if not findings:
        return None
    triads = [row for row in findings if row["type"] == "catalytic_triad"]
    pockets = [row for row in findings if row["type"] == "pocket_ph_sensitive"]
    parts: list[str] = []
    if triads:
        parts.append(
            "catalytic triad(s): " + "; ".join(
                "-".join(row["residues"]) for row in triads))
    if pockets:
        names = ", ".join(
            residue for row in pockets for residue in row["residues"])
        parts.append(f"pH-sensitive residues near the pocket: {names}")
    return (
        "Catalytic residues detected near the binding site; protonation "
        "states may be wrong without pKa analysis. Hydrogens were added "
        "by toolkit template (pH 7.4 nominal), NOT by a pKa calculation "
        f"({'; '.join(parts)}). Run PROPKA / H++ and pass the "
        "protonated PDB as the receptor input for rigorous work."
    )
