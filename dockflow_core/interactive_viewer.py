"""Interactive 3D viewer HTML generation (peer item 20).

Writes ``<run>/visualization/interactive.html``: a single self-contained
page (3Dmol.js from CDN, no Python-side install) that renders

* the prepared receptor (cartoon + stick hetero atoms),
* every top pose of every ligand as sticks, with a pose selector,
* the grid box (wireframe),
* detected contact residues as clickable spheres with labels,
* the co-crystallized reference ligand overlaid as semi-transparent
  sticks when redocking validation data exists.

Everything the page needs is inlined as one JSON payload, so the file
can be archived with the run and opened offline (the only network use
is the 3Dmol.js CDN script; users behind a firewall can swap the script
tag for a local copy - noted in the page footer).
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from .utils import get_logger

logger = get_logger("viewer3d")

_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>DockFlow interactive viewer - __RUN_ID__</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
 body { margin: 0; font-family: system-ui, sans-serif; background: #f7f8fa; }
 #layout { display: flex; height: 100vh; }
 #viewport { flex: 1; position: relative; }
 #panel { width: 300px; padding: 12px 16px; overflow-y: auto;
          border-left: 1px solid #ddd; background: #fff; }
 h1 { font-size: 15px; margin: 4px 0 10px; }
 .muted { color: #666; font-size: 12px; }
 button { margin: 2px 2px 6px 0; padding: 4px 8px; font-size: 12px;
          cursor: pointer; border: 1px solid #bbb; border-radius: 4px;
          background: #fff; }
 button.active { background: #1a7f37; color: #fff; border-color: #1a7f37; }
 #info { font-size: 12px; border-top: 1px solid #ddd; padding-top: 8px;
         margin-top: 8px; min-height: 60px; }
 table { border-collapse: collapse; font-size: 11px; }
 td, th { border: 1px solid #eee; padding: 2px 5px; text-align: left; }
</style>
</head>
<body>
<div id="layout">
  <div id="viewport"></div>
  <div id="panel">
    <h1>__RUN_ID__</h1>
    <div class="muted">receptor: __ENGINE__ | backend: __BACKEND__</div>
    <h3>poses</h3><div id="poses"></div>
    <h3>grid box</h3>
    <div class="muted">center __CENTER__ A, size __SIZE__ A
      (source: __BOX_SOURCE__)</div>
    <label><input type="checkbox" id="box" checked> show box</label>
    <label><input type="checkbox" id="ref" checked> crystal reference</label>
    <h3>contact residues (top pose)</h3>
    <table id="contacts"></table>
    <div id="info" class="muted">click a contact residue for details</div>
    <div class="muted" style="margin-top:10px">
      heuristic geometric contacts (distance + atom type only);
      see docs/interaction_criteria.md. Offline: replace the 3Dmol.js
      script tag with a local copy.
    </div>
  </div>
</div>
<script>
const DATA = __PAYLOAD__;
const viewer = $3Dmol.createViewer("viewport",
    {backgroundColor: "#f7f8fa", antialias: true});

// receptor: model 0
viewer.addModel(DATA.receptor, "pdbqt");
viewer.setStyle({model: 0}, {cartoon: {color: "spectrum"}});
viewer.setStyle({model: 0, atom: ["CA", "N", "O", "C"]}, {});
viewer.addStyle({model: 0, hetflag: true},
                {stick: {radius: 0.15, colorscheme: "greenCarbon"}});

// ligand poses (one model each)
DATA.poses.forEach((pose, index) => {
  viewer.addModel(pose.pdbqt, "pdbqt");
  viewer.setStyle({model: index + 1},
                  {stick: {radius: 0.16, colorscheme: "cyanCarbon"},
                   sphere: {scale: 0.14, opacity: 0.35}});
});

// crystal reference overlay (semi-transparent sticks)
if (DATA.reference) {
  viewer.addModel(DATA.reference, "pdb");
  viewer.setStyle({model: DATA.poses.length + 1},
                  {stick: {radius: 0.14, opacity: 0.45,
                           colorscheme: "magentaCarbon"}});
}

// grid box wireframe
function drawBox(center, size) {
  const h = [size[0] / 2, size[1] / 2, size[2] / 2];
  const c = center;
  const corners = [];
  for (const sx of [-1, 1]) for (const sy of [-1, 1])
    for (const sz of [-1, 1])
      corners.push([c[0] + sx * h[0], c[1] + sy * h[1], c[2] + sz * h[2]]);
  const edges = [[0,1],[0,2],[0,4],[1,3],[1,5],[2,3],[2,6],[3,7],
                 [4,5],[4,6],[5,7],[6,7]];
  edges.forEach(([a, b]) => viewer.addLine(
      {start: {x: corners[a][0], y: corners[a][1], z: corners[a][2]},
       end:   {x: corners[b][0], y: corners[b][1], z: corners[b][2]},
       color: "#c62828"}));
}
drawBox(DATA.box.center, DATA.box.size);

// contact residues as clickable spheres
const contactAtoms = {};
DATA.contacts.forEach((contact, index) => {
  viewer.addSphere({center: {x: contact.xyz[0], y: contact.xyz[1],
                             z: contact.xyz[2]},
                    radius: 1.4, color: "#f5b301", alpha: 0.35});
  contactAtoms[index] = contact;
});
viewer.addEventListener("click",
    (event, viewerImpl) => { viewer.render(); });

// pose selector UI
const poseDiv = document.getElementById("poses");
let activePose = 0;
function showPose(index) {
  activePose = index;
  DATA.poses.forEach((pose, i) => {
    viewer.getModel(i + 1).show(i === index);
  });
  document.querySelectorAll("#poses button").forEach(
      (button, i) => button.classList.toggle("active", i === index));
  const pose = DATA.poses[index];
  document.getElementById("info").innerHTML =
      "<b>" + pose.ligand + " pose " + pose.pose + "</b><br>" +
      "affinity " + pose.affinity + " kcal/mol<br>" +
      (pose.crystal_rmsd != null
        ? "crystal RMSD " + pose.crystal_rmsd + " A (" +
          DATA.rmsd_method + ")" : "no reference RMSD");
  viewer.render();
}
DATA.poses.forEach((pose, index) => {
  const button = document.createElement("button");
  button.textContent = pose.ligand + " #" + pose.pose +
      " (" + pose.affinity + ")";
  button.onclick = () => showPose(index);
  poseDiv.appendChild(button);
});

// contact table (clickable rows)
const contactTable = document.getElementById("contacts");
DATA.contacts.forEach((contact, index) => {
  const row = contactTable.insertRow(-1);
  const cell = row.insertCell(0);
  cell.textContent = contact.resname + contact.resseq +
      " (" + (contact.chain || "-") + ")";
  const cell2 = row.insertCell(1);
  cell2.textContent = contact.kinds;
  const cell3 = row.insertCell(2);
  cell3.textContent = contact.closest;
  row.onclick = () => {
    document.getElementById("info").innerHTML =
        "<b>" + contact.resname + contact.resseq + " chain " +
        contact.chain + "</b><br>" + contact.kinds +
        "<br>closest ligand atom " + contact.closest + " A";
    viewer.zoomTo({resi: contact.resseq, chain: contact.chain});
    viewer.render();
  };
});

document.getElementById("box").onchange = (event) => {
  viewer.render();
};
document.getElementById("ref").onchange = (event) => {
  const refModel = viewer.getModel(DATA.poses.length + 1);
  if (refModel) { refModel.show(event.target.checked); viewer.render(); }
};

showPose(0);
viewer.zoomTo({model: 1});
viewer.render();
</script>
</body>
</html>
"""


def _pdbqt_text(path: Path | None) -> str:
    if path is None or not path.is_file():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def collect_pose_files(run_dir: Path, ligand_name: str) -> list[Path]:
    """Pose PDBQT files of one ligand (split per-model files, else the
    multi-model output)."""
    analysis_poses = sorted(
        (run_dir / "analysis").glob(f"{ligand_name}_pose*.pdbqt"))
    if analysis_poses:
        return analysis_poses
    combined = run_dir / "docking" / f"{ligand_name}_out.pdbqt"
    return [combined] if combined.is_file() else []


def build_interactive_html(run_dir: str | Path) -> Path | None:
    """Generate ``visualization/interactive.html`` for a completed run.

    Returns the path (or ``None`` when the run has no docking outputs
    yet - e.g. a failed run - so the viewer is simply skipped).
    """
    run_dir = Path(run_dir)
    docking_dir = run_dir / "docking"
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        logger.debug("no manifest; interactive viewer skipped")
        return None
    manifest: dict[str, Any] = json.loads(
        manifest_path.read_text(encoding="utf-8"))
    results = (manifest.get("docking") or {}).get("results") or []
    ok_results = [entry for entry in results if entry.get("poses")]

    receptor = Path((manifest.get("receptor") or {}).get("pdbqt", ""))
    receptor_text = _pdbqt_text(receptor if receptor.is_absolute() else
                                run_dir / "prepared" / "receptor.pdbqt")
    if not receptor_text:
        logger.debug("no receptor PDBQT; interactive viewer skipped")
        return None

    poses: list[dict[str, Any]] = []
    for entry in ok_results:
        ligand = entry.get("ligand_name", "ligand")
        for pose in entry.get("poses", [])[:5]:
            poses.append({
                "ligand": ligand,
                "pose": pose.get("model", 1),
                "affinity": pose.get("affinity"),
                "crystal_rmsd": pose.get("crystal_rmsd"),
                "pdbqt": _split_pose_model(
                    _pdbqt_text(docking_dir / f"{ligand}_out.pdbqt"),
                    int(pose.get("model", 1))),
            })
    if not poses:
        logger.debug("no poses; interactive viewer skipped")
        return None

    gridbox = manifest.get("gridbox") or {}
    box = {
        "center": [round(float(v), 2) for v in gridbox.get("center", (0, 0, 0))],
        "size": [round(float(v), 2) for v in gridbox.get("size", (0, 0, 0))],
    }

    # contacts of the best ligand's top pose
    analysis_payload = ((manifest.get("analysis") or {}).get("poses") or {})
    best_ligand = poses[0]["ligand"]
    contacts = []
    for pose_entry in analysis_payload.get(best_ligand, [])[:1]:
        for residue in pose_entry.get("residues", [])[:40]:
            contacts.append({
                "chain": residue.get("chain", ""),
                "resname": residue.get("resname", ""),
                "resseq": residue.get("resseq", 0),
                "kinds": (f"{residue.get('geom_hbond_contacts', 0)} hbond, "
                          f"{residue.get('hydrophobic_contacts', 0)} hydro, "
                          f"{residue.get('ionic_contacts', 0)} ionic, "
                          f"{residue.get('metal_contacts', 0)} metal"),
                "closest": residue.get("closest"),
                "xyz": None,
            })
    # approximate contact sphere positions from the receptor CA atoms
    if contacts:
        from .pdbio import parse_pdbqt

        try:
            ca_atoms = {
                (a.chain.strip(), int(a.resseq)): (a.x, a.y, a.z)
                for a in parse_pdbqt(receptor).atoms
                if (a.name or "").strip().upper() == "CA"
            }
            for contact in contacts:
                key = (contact["chain"], int(contact["resseq"]))
                contact["xyz"] = list(ca_atoms.get(key, (0, 0, 0)))
        except Exception:  # noqa: BLE001 - spheres are decorative
            for contact in contacts:
                contact["xyz"] = [0, 0, 0]

    # crystal reference overlay (item 20c)
    reference_text = ""
    reference = (manifest.get("paths") or {}).get("crystal_rmsd_reference")
    if reference and " in " in str(reference):
        resname, structure = str(reference).split(" in ", 1)
        structure_path = run_dir / "raw" / structure
        if structure_path.is_file():
            try:
                from .analyzer import extract_reference_atoms

                atoms = extract_reference_atoms(structure_path,
                                                resname.strip())
                reference_text = _atoms_to_pdb(atoms)
            except Exception:  # noqa: BLE001 - overlay is optional
                logger.debug("reference overlay unavailable")

    payload = {
        "run_id": manifest.get("run_id", run_dir.name),
        "receptor": receptor_text,
        "poses": [pose for pose in poses if pose["pdbqt"]],
        "box": box,
        "box_source": gridbox.get("source", "?"),
        "contacts": [c for c in contacts if c["xyz"]],
        "reference": reference_text,
        "rmsd_method": (manifest.get("analysis") or {}).get(
            "crystal_rmsd_method") or "symmetry-aware RMSD",
    }
    if not payload["poses"]:
        logger.debug("no pose PDBQT data; interactive viewer skipped")
        return None

    page = _TEMPLATE
    replacements = {
        "__RUN_ID__": html.escape(str(manifest.get("run_id", run_dir.name))),
        "__ENGINE__": html.escape(
            str((manifest.get("receptor") or {}).get("engine", "?"))),
        "__BACKEND__": html.escape(
            str((manifest.get("docking") or {}).get("backend", "?"))),
        "__CENTER__": str(box["center"]),
        "__SIZE__": str(box["size"]),
        "__BOX_SOURCE__": html.escape(str(gridbox.get("source", "?"))),
        "__PAYLOAD__": json.dumps(payload),
    }
    for key, value in replacements.items():
        page = page.replace(key, value)
    out_path = run_dir / "visualization" / "interactive.html"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(page, encoding="utf-8")
    logger.info("interactive viewer: %s", out_path)
    return out_path


def _split_pose_model(pdbqt_text: str, model_index: int) -> str:
    """Extract one MODEL block from a multi-model PDBQT string."""
    if not pdbqt_text:
        return ""
    lines = pdbqt_text.splitlines()
    current = 0
    selected: list[str] = []
    in_model = False
    for line in lines:
        if line.startswith("MODEL"):
            current += 1
            in_model = current == model_index
            continue
        if line.startswith("ENDMDL"):
            if in_model:
                break
            in_model = False
            continue
        if in_model and line.startswith(("ATOM", "HETATM")):
            selected.append(line)
    if not selected and model_index == 1:
        # single-model file without MODEL blocks
        selected = [line for line in lines
                    if line.startswith(("ATOM", "HETATM"))]
    return "\n".join(selected) + "\n" if selected else ""


def _atoms_to_pdb(atoms) -> str:
    """Minimal PDB block for the reference overlay atoms."""
    lines = []
    for index, atom in enumerate(atoms, start=1):
        element = (atom.element or "C").strip()[:2]
        name = (atom.name or "X").strip()[:4]
        lines.append(
            f"HETATM{index:5d} {name:<4s} LIG A   1    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}  "
        )
    lines.append("END")
    return "\n".join(lines) + "\n"
