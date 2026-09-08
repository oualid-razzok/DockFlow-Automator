#!/usr/bin/env python3
"""Fix conda-build placeholder prefixes in an extracted MGLTools tree.

The bioconda ``mgltools 1.5.7`` package (the last release, 2011,
Python-2-based) embeds its build-time placeholder prefix in scripts
(shebangs / ``MGL_ROOT``) and ELF binaries (RPATH / compiled-in
prefixes).  Extracting the tarball does not relocate them, so
``prepare_receptor4.py`` fails with ``not found`` errors until the
prefixes are rewritten.

Text files get a plain replacement; binaries get a NUL-padded,
length-preserving replacement (byte-safe).  Run this ONCE after
extracting the package, then point ``compare_workflow.py --mgltools``
at the tree.

Usage::

    python scripts/fix_mgltools_prefix.py --root /path/to/mgltools
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True,
                        help="extracted mgltools tree (contains bin/)")
    args = parser.parse_args(argv)

    root = Path(args.root).resolve()
    if not (root / "bin" / "prepare_receptor4.py").is_file():
        print(f"error: {root} does not look like an MGLTools tree "
              "(bin/prepare_receptor4.py missing)", file=sys.stderr)
        return 2

    placeholder_re = re.compile(
        rb"/opt/conda/conda-bld/mgltools_[0-9]+/_h_env_placehold[a-z_]*")
    new_prefix = str(root).encode("utf-8")

    probe = (root / "bin" / "pythonsh").read_bytes()
    match = placeholder_re.search(probe)
    if match is None:
        print("no conda placeholder prefix found - tree already relocated")
        return 0
    placeholder = match.group(0)
    if len(placeholder) < len(new_prefix):
        print(f"error: placeholder ({len(placeholder)} chars) is shorter "
              f"than the new prefix ({len(new_prefix)} chars); extract the "
              "tree to a shorter path", file=sys.stderr)
        return 3

    def is_elf(path: Path) -> bool:
        try:
            with open(path, "rb") as fh:
                return fh.read(4) == b"\x7fELF"
        except OSError:
            return False

    def is_text(path: Path) -> bool:
        try:
            with open(path, "rb") as fh:
                return b"\x00" not in fh.read(4096)
        except OSError:
            return False

    patched_text = patched_bin = 0
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix in (".tar", ".bz2", ".conda"):
            continue
        try:
            data = path.read_bytes()
        except OSError:
            continue
        if placeholder not in data:
            continue
        if is_elf(path) or not is_text(path):
            replacement = new_prefix + b"\x00" * (
                len(placeholder) - len(new_prefix))
            path.write_bytes(data.replace(placeholder, replacement))
            patched_bin += 1
        else:
            path.write_bytes(data.replace(placeholder, new_prefix))
            patched_text += 1
    print(f"relocated {patched_text} script(s) and {patched_bin} "
          f"binary file(s) to {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
