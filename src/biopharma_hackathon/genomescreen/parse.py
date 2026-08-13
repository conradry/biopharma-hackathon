"""Parsers for GenomeScreen's filename and grid-file conventions.

Every screen result lives in a directory whose name encodes the target::

    AF-<uniprot_acc>-F1-model_v<version>_<fragment_idx>_<pocket>

``<pocket>`` is either ``pocket<N>`` (a pocket found by apo pocket detection) or a bare
``<N>`` (a pocket transferred from an aligned holo template).  Inside, each candidate
receptor conformation is a pair of files::

    <pocket_key>_<structure_idx>[_<template_pdb_id>]_complex_refined.pdbgz
    <pocket_key>_<structure_idx>[_<template_pdb_id>]_complex_refined_grid.in

The ``<template_pdb_id>`` segment is present only for template-derived pockets, where it
records the PDB entry the pocket was copied from.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

#: Directory name -> UniProt accession, AlphaFold model version, fragment, pocket.
POCKET_DIR_RE = re.compile(
    r"^AF-(?P<acc>[A-Z0-9]+)-F(?P<entry>\d+)-model_v(?P<model_version>\d+)"
    r"_(?P<fragment_idx>\d+)_(?P<pocket>pocket(?P<detected_idx>\d+)|(?P<template_idx>\d+))$"
)

#: The part of a structure filename that follows the pocket key.
STRUCTURE_SUFFIX_RE = re.compile(
    r"^(?P<structure_idx>\d+)(?:_(?P<template_pdb_id>[0-9a-z][0-9a-z]{3}))?_complex_refined$"
)

STRUCTURE_EXT = ".pdbgz"
GRID_EXT = "_grid.in"
HITS_FILENAME = "leader.csv"

POCKET_KIND_DETECTED = "detected"
POCKET_KIND_TEMPLATE = "template"


@dataclass(frozen=True, slots=True)
class PocketId:
    """A parsed screen-result directory name."""

    pocket_key: str
    uniprot_acc: str
    af_entry: int
    af_model_version: int
    fragment_idx: int
    pocket_kind: str
    pocket_idx: int


@dataclass(frozen=True, slots=True)
class StructureFiles:
    """One refined receptor conformation and its docking grid."""

    pocket_key: str
    structure_idx: int
    template_pdb_id: str | None
    structure_path: str
    grid_path: str | None
    grid_center_x: float | None
    grid_center_y: float | None
    grid_center_z: float | None
    grid_file: str | None
    receptor_file: str | None


def parse_pocket_dirname(dirname: str) -> PocketId:
    """Parse a screen-result directory name.

    Raises:
        ValueError: if the name does not follow the GenomeScreen convention.
    """
    match = POCKET_DIR_RE.match(dirname)
    if match is None:
        raise ValueError(f"not a GenomeScreen pocket directory name: {dirname!r}")
    detected = match.group("detected_idx")
    kind = POCKET_KIND_DETECTED if detected is not None else POCKET_KIND_TEMPLATE
    return PocketId(
        pocket_key=dirname,
        uniprot_acc=match.group("acc"),
        af_entry=int(match.group("entry")),
        af_model_version=int(match.group("model_version")),
        fragment_idx=int(match.group("fragment_idx")),
        pocket_kind=kind,
        pocket_idx=int(detected if detected is not None else match.group("template_idx")),
    )


def parse_grid_file(path: str | os.PathLike[str]) -> dict[str, object]:
    """Read a Schrodinger-style ``*_grid.in`` file.

    Returns a dict with ``grid_file``, ``receptor_file`` and ``grid_center`` (a 3-tuple of
    floats); any key absent from the file is omitted.
    """
    out: dict[str, object] = {}
    with open(path) as handle:
        for line in handle:
            key, _, rest = line.strip().partition(" ")
            rest = rest.strip()
            if not rest:
                continue
            if key == "GRIDFILE":
                out["grid_file"] = rest
            elif key == "RECEP_FILE":
                out["receptor_file"] = rest
            elif key == "GRID_CENTER":
                x, y, z = (float(part) for part in rest.split(","))
                out["grid_center"] = (x, y, z)
    return out


def scan_pocket_dir(directory: str | os.PathLike[str]) -> tuple[PocketId, list[StructureFiles]]:
    """Parse one screen-result directory into a pocket and its receptor structures.

    Structures whose ``*_grid.in`` file is missing (three exist in the public release) are
    still returned, with all grid fields set to ``None``.
    """
    directory = Path(directory)
    pocket = parse_pocket_dirname(directory.name)
    prefix = f"{pocket.pocket_key}_"

    structures: list[str] = []
    grids: set[str] = set()
    for entry in os.scandir(directory):
        name = entry.name
        if name == HITS_FILENAME:
            continue
        if name.endswith(GRID_EXT):
            grids.add(name[: -len(GRID_EXT)])
        elif name.endswith(STRUCTURE_EXT):
            structures.append(name[: -len(STRUCTURE_EXT)])
        else:
            raise ValueError(f"unexpected file in {directory}: {name!r}")

    parsed: list[StructureFiles] = []
    for stem in structures:
        if not stem.startswith(prefix):
            raise ValueError(f"structure {stem!r} does not belong to pocket {pocket.pocket_key!r}")
        match = STRUCTURE_SUFFIX_RE.match(stem[len(prefix) :])
        if match is None:
            raise ValueError(f"unparseable structure filename: {stem!r}")

        grid_path = directory / f"{stem}{GRID_EXT}" if stem in grids else None
        grid = parse_grid_file(grid_path) if grid_path is not None else {}
        center = grid.get("grid_center")
        cx, cy, cz = center if isinstance(center, tuple) else (None, None, None)
        parsed.append(
            StructureFiles(
                pocket_key=pocket.pocket_key,
                structure_idx=int(match.group("structure_idx")),
                template_pdb_id=match.group("template_pdb_id"),
                structure_path=str(directory / f"{stem}{STRUCTURE_EXT}"),
                grid_path=str(grid_path) if grid_path is not None else None,
                grid_center_x=cx,
                grid_center_y=cy,
                grid_center_z=cz,
                grid_file=grid.get("grid_file"),  # type: ignore[arg-type]
                receptor_file=grid.get("receptor_file"),  # type: ignore[arg-type]
            )
        )

    parsed.sort(key=lambda s: s.structure_idx)
    return pocket, parsed
