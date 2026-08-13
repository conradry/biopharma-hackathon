"""Ingest the GenomeScreen virtual-screening dataset into DuckDB.

GenomeScreen (Jia et al. 2024, https://doi.org/10.1101/2024.09.02.610777) ships as a
directory of per-pocket screen results.  This package turns that tree into a single
queryable DuckDB file; see :func:`build_database`.
"""

from biopharma_hackathon.genomescreen.ingest import build_database, connect
from biopharma_hackathon.genomescreen.parse import (
    PocketId,
    StructureFiles,
    parse_grid_file,
    parse_pocket_dirname,
    scan_pocket_dir,
)

__all__ = [
    "PocketId",
    "StructureFiles",
    "build_database",
    "connect",
    "parse_grid_file",
    "parse_pocket_dirname",
    "scan_pocket_dir",
]
