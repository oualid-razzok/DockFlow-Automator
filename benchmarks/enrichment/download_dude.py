"""Download + prepare a DUD-E actives/decoys panel (audit item 22).

Fetches SMILES lists for one DUD-E target, prepares every molecule into
PDBQT (DockFlow ligand preparation), and writes an actives id list for
`dockflow enrich`.

Usage::

    python benchmarks/enrichment/download_dude.py --target hivpr --out data
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

DUDE_BASE = "https://dude.docking.org/targets/{target}"
FILES = {"actives": "actives_final.ism", "decoys": "decoys_final.ism"}


def fetch(target: str, out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for kind, filename in FILES.items():
        url = f"{DUDE_BASE.format(target=target)}/{filename}"
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        path = out_dir / filename
        path.write_text(response.text, encoding="utf-8")
        paths[kind] = path
        print(f"downloaded {len(response.text.splitlines())} {kind} -> {path}")
    return paths


def prepare_panel(paths: dict[str, Path], out_dir: Path) -> Path:
    """Merge actives+decoys into one multi-SMILES file and prepare it."""
    merged = out_dir / "panel.smi"
    actives_ids: list[str] = []
    lines: list[str] = []
    for kind in ("actives", "decoys"):
        for line in paths[kind].read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) < 2:
                continue
            smiles, molecule_id = parts[0], parts[1]
            lines.append(f"{smiles} {molecule_id}")
            if kind == "actives":
                actives_ids.append(molecule_id)
    merged.write_text("\n".join(lines) + "\n", encoding="utf-8")
    actives_list = out_dir / "actives.txt"
    actives_list.write_text("\n".join(actives_ids) + "\n", encoding="utf-8")
    print(f"panel: {len(lines)} molecules ({len(actives_ids)} actives)")
    print(f"next: dockflow prep ligand --in {merged} --library "
          f"--out-dir {out_dir / 'prepared'}")
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="hivpr",
                        help="DUD-E target directory name")
    parser.add_argument("--out", default="data")
    args = parser.parse_args()
    out_dir = Path(args.out)
    paths = fetch(args.target, out_dir)
    prepare_panel(paths, out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
