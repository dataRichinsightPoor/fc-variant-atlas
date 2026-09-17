#!/usr/bin/env python3
"""Build a standalone installable zip of the PyMOL plugin.

The plugin reads the atlas from the first data directory it finds, and an
installed copy sits outside the repository, so the data files are copied into
the package before zipping. That makes the zip self-contained: a reader with
PyMOL and no Python packaging knowledge can install it from the Plugin Manager
and get the same numbers this repository publishes.

    python pymol_plugin/build_plugin_zip.py

writes dist/fcatlas_pymol-<version>.zip
"""

from __future__ import annotations

import re
import shutil
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PACKAGE = HERE / "fcatlas_pymol"
DATA = REPO / "data"
DIST = REPO / "dist"

#: What the plugin needs at runtime, and nothing else.
BUNDLED = ("atlas.json", "interface_detail.json", "structures.json")


def version() -> str:
    text = (PACKAGE / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'__version__ = "([^"]+)"', text)
    if not match:
        raise SystemExit("no __version__ in the plugin")
    return match.group(1)


def stage_data() -> list[Path]:
    target = PACKAGE / "data"
    target.mkdir(exist_ok=True)
    staged = []
    for name in BUNDLED:
        source = DATA / name
        if not source.is_file():
            raise SystemExit(f"missing {source}; regenerate the atlas first")
        destination = target / name
        shutil.copy2(source, destination)
        staged.append(destination)
    return staged


def build() -> Path:
    staged = stage_data()
    DIST.mkdir(exist_ok=True)
    out = DIST / f"fcatlas_pymol-{version()}.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(PACKAGE.rglob("*")):
            if path.is_dir() or "__pycache__" in path.parts:
                continue
            archive.write(path, path.relative_to(HERE))
        archive.write(HERE / "README.md", "fcatlas_pymol/README.md")
        archive.write(REPO / "LICENSE", "fcatlas_pymol/LICENSE")
    sizes = ", ".join(f"{p.name} {p.stat().st_size // 1024} kB" for p in staged)
    print(f"{out.relative_to(REPO)}  ({out.stat().st_size // 1024} kB)")
    print(f"bundled data: {sizes}")
    return out


def verify(path: Path) -> None:
    """Check the zip has what the plugin manager and the plugin both need."""
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
    required = {"fcatlas_pymol/__init__.py"} | {
        f"fcatlas_pymol/data/{name}" for name in BUNDLED
    }
    missing = required - names
    if missing:
        raise SystemExit("incomplete zip, missing: " + ", ".join(sorted(missing)))
    print(f"verified {len(names)} entries")


if __name__ == "__main__":
    archive = build()
    verify(archive)
    if "--keep-staged" not in sys.argv:
        shutil.rmtree(PACKAGE / "data", ignore_errors=True)
        print("staged data removed from the working tree")
