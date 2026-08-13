"""Build a queryable DuckDB database from a GenomeScreen release directory.

Usage::

    uv run genomescreen-ingest /home/ubuntu/datasets/GenomeScreenDB genomescreen.duckdb

The 27 GB of ``.pdbgz`` receptor structures stay on disk; the database records their paths
in ``pocket_structures``.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import time
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb

from biopharma_hackathon.genomescreen.parse import HITS_FILENAME, scan_pocket_dir

SCHEMA_PATH = Path(__file__).with_name("schema.sql")

#: Threads used to walk pocket directories.  The scan is I/O bound on ~215k small grid
#: files, so oversubscribing cores pays off.
DEFAULT_SCAN_WORKERS = min(32, (os.cpu_count() or 4) * 4)

#: Directory inside a release that holds one subdirectory per screened pocket.
SCREEN_RESULTS_DIRNAME = "screen_results"

#: Molecule id namespaces, tried in order.  Regexes are DuckDB (RE2) syntax.
SOURCE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("ZINC", r"ZINC\d{12}"),
    ("Enamine REAL", r"Z\d+(_\d+)?(_T\d+)?"),
    ("Enamine PV", r"PV-\d+(_\d+)?(_T\d+)?"),
)


def resolve_screen_results(root: str | os.PathLike[str]) -> Path:
    """Accept either a release root or its ``screen_results`` directory."""
    root = Path(root)
    if root.name == SCREEN_RESULTS_DIRNAME:
        return root
    nested = root / SCREEN_RESULTS_DIRNAME
    if nested.is_dir():
        return nested
    raise FileNotFoundError(f"no {SCREEN_RESULTS_DIRNAME}/ under {root}")


def iter_pocket_dirs(screen_results: Path, limit: int | None = None) -> Iterator[Path]:
    """Yield pocket directories in sorted order, for reproducible builds."""
    names = sorted(entry.name for entry in os.scandir(screen_results) if entry.is_dir())
    if limit is not None:
        names = names[:limit]
    for name in names:
        yield screen_results / name


def _source_case_expr(column: str) -> str:
    branches = "\n        ".join(
        f"WHEN regexp_full_match({column}, '{pattern}') THEN '{name}'"
        for name, pattern in SOURCE_PATTERNS
    )
    return f"CASE\n        {branches}\n        ELSE 'unknown'\n    END"


_POCKET_COLUMNS = {
    "pocket_key": "VARCHAR",
    "uniprot_acc": "VARCHAR",
    "af_entry": "INTEGER",
    "af_model_version": "INTEGER",
    "fragment_idx": "INTEGER",
    "pocket_kind": "VARCHAR",
    "pocket_idx": "INTEGER",
    "dir_path": "VARCHAR",
}

_STRUCTURE_COLUMNS = {
    "pocket_key": "VARCHAR",
    "structure_idx": "INTEGER",
    "template_pdb_id": "VARCHAR",
    "structure_path": "VARCHAR",
    "grid_path": "VARCHAR",
    "grid_center_x": "DOUBLE",
    "grid_center_y": "DOUBLE",
    "grid_center_z": "DOUBLE",
    "grid_file": "VARCHAR",
    "receptor_file": "VARCHAR",
}


def _copy_csv_into(
    conn: duckdb.DuckDBPyConnection, table: str, path: Path, columns: dict[str, str]
) -> None:
    """Bulk-load a headerless staging CSV.

    DuckDB's per-statement overhead makes ``executemany`` of a few hundred thousand rows
    take minutes; going through its CSV reader is two orders of magnitude faster.
    """
    spec = ", ".join(f"'{name}': '{sql_type}'" for name, sql_type in columns.items())
    conn.execute(
        f"INSERT INTO {table} ({', '.join(columns)}) "
        "SELECT * FROM read_csv(?, header = false, columns = {" + spec + "})",
        [str(path)],
    )


def _load_pocket_metadata(
    conn: duckdb.DuckDBPyConnection,
    screen_results: Path,
    limit: int | None,
    workers: int,
) -> tuple[int, int]:
    """Populate ``pockets`` and ``pocket_structures`` by walking the directory tree."""
    directories = list(iter_pocket_dirs(screen_results, limit))
    n_pockets = 0
    n_structures = 0

    with tempfile.TemporaryDirectory(prefix="genomescreen-stage-") as staging:
        pockets_csv = Path(staging) / "pockets.csv"
        structures_csv = Path(staging) / "structures.csv"
        with (
            open(pockets_csv, "w", newline="") as pf,
            open(structures_csv, "w", newline="") as sf,
            ThreadPoolExecutor(max_workers=workers) as pool,
        ):
            pocket_writer = csv.writer(pf)
            structure_writer = csv.writer(sf)
            # map() preserves input order, so a rebuild of the same tree is byte-identical.
            for pocket, structures in pool.map(scan_pocket_dir, directories):
                pocket_writer.writerow(
                    (
                        pocket.pocket_key,
                        pocket.uniprot_acc,
                        pocket.af_entry,
                        pocket.af_model_version,
                        pocket.fragment_idx,
                        pocket.pocket_kind,
                        pocket.pocket_idx,
                        str(screen_results / pocket.pocket_key),
                    )
                )
                structure_writer.writerows(
                    (
                        s.pocket_key,
                        s.structure_idx,
                        s.template_pdb_id,
                        s.structure_path,
                        s.grid_path,
                        s.grid_center_x,
                        s.grid_center_y,
                        s.grid_center_z,
                        s.grid_file,
                        s.receptor_file,
                    )
                    for s in structures
                )
                n_pockets += 1
                n_structures += len(structures)

        if n_pockets:
            _copy_csv_into(conn, "pockets", pockets_csv, _POCKET_COLUMNS)
        if n_structures:
            _copy_csv_into(conn, "pocket_structures", structures_csv, _STRUCTURE_COLUMNS)

    return n_pockets, n_structures


def _load_hits(conn: duckdb.DuckDBPyConnection, screen_results: Path) -> None:
    """Read every ``leader.csv`` at once and split it into ``molecules`` and ``hits``.

    Only pockets already present in ``pockets`` are kept, so ``--limit`` builds stay
    self-consistent.
    """
    glob = str(screen_results / "*" / HITS_FILENAME)
    conn.execute(
        """
        CREATE TEMP TABLE raw_hits AS
        SELECT
            regexp_extract(filename, '([^/]+)/[^/]+$', 1) AS pocket_key,
            smiles,
            oid,
            score,
            "Name" AS source_index
        FROM read_csv(
            ?,
            header = true,
            columns = {
                'smiles': 'VARCHAR',
                'Name': 'BIGINT',
                'oid': 'VARCHAR',
                'score': 'DOUBLE'
            },
            filename = true
        )
        """,
        [glob],
    )
    conn.execute("DELETE FROM raw_hits WHERE pocket_key NOT IN (SELECT pocket_key FROM pockets)")

    # An oid can appear with more than one SMILES (protomer/tautomer variants of the same
    # compound), so the molecule grain is the (oid, smiles) pair.  catalog_id collapses them.
    conn.execute(f"""
        INSERT INTO molecules
            (mol_id, oid, smiles, source, catalog_id, protomer_idx, tautomer_idx)
        SELECT
            row_number() OVER (ORDER BY oid, smiles) - 1 AS mol_id,
            oid,
            smiles,
            source,
            CASE WHEN source = 'ZINC' THEN oid ELSE regexp_extract(oid, '^([^_]+)', 1) END,
            try_cast(regexp_extract(oid, '^[^_]+_(\\d+)', 1) AS INTEGER),
            try_cast(regexp_extract(oid, '_T(\\d+)$', 1) AS INTEGER)
        FROM (
            SELECT DISTINCT oid, smiles, {_source_case_expr("oid")} AS source FROM raw_hits
        )
    """)

    conn.execute("""
        INSERT INTO hits (pocket_key, mol_id, score, rank_in_pocket, source_index)
        SELECT
            r.pocket_key,
            m.mol_id,
            r.score,
            row_number() OVER (
                PARTITION BY r.pocket_key ORDER BY r.score DESC, m.mol_id
            ) AS rank_in_pocket,
            r.source_index
        FROM raw_hits AS r
        JOIN molecules AS m ON m.oid = r.oid AND m.smiles = r.smiles
    """)
    conn.execute("DROP TABLE raw_hits")


def _fill_rollups(conn: duckdb.DuckDBPyConnection) -> None:
    conn.execute("""
        UPDATE pockets AS p
        SET n_hits = agg.n_hits, best_score = agg.best_score
        FROM (
            SELECT pocket_key, count(*) AS n_hits, max(score) AS best_score
            FROM hits GROUP BY 1
        ) AS agg
        WHERE agg.pocket_key = p.pocket_key
    """)
    conn.execute("""
        UPDATE pockets AS p
        SET n_structures = agg.n
        FROM (SELECT pocket_key, count(*) AS n FROM pocket_structures GROUP BY 1) AS agg
        WHERE agg.pocket_key = p.pocket_key
    """)
    conn.execute("""
        INSERT INTO targets (uniprot_acc, n_pockets, n_hits, best_score)
        SELECT uniprot_acc, count(*), sum(n_hits), max(best_score)
        FROM pockets GROUP BY 1
    """)


def build_database(
    root: str | os.PathLike[str],
    database: str | os.PathLike[str],
    *,
    limit: int | None = None,
    overwrite: bool = False,
    threads: int | None = None,
    scan_workers: int = DEFAULT_SCAN_WORKERS,
    progress: bool = False,
) -> Path:
    """Build ``database`` from the GenomeScreen release at ``root``.

    Args:
        root: release directory (or its ``screen_results`` subdirectory).
        database: path of the DuckDB file to create.
        limit: ingest only the first N pocket directories, for smoke tests.
        overwrite: replace ``database`` if it already exists.
        threads: DuckDB thread count; ``None`` uses DuckDB's default.
        scan_workers: threads used for the directory walk.
        progress: print per-stage timings to stderr.

    Returns:
        The path to the database that was written.
    """
    screen_results = resolve_screen_results(root)
    database = Path(database)
    if database.exists():
        if not overwrite:
            raise FileExistsError(f"{database} already exists; pass overwrite=True to replace it")
        database.unlink()

    def stage(message: str, started: float) -> None:
        if progress:
            print(f"  {message} ({time.monotonic() - started:.1f}s)", file=sys.stderr)

    conn = duckdb.connect(str(database))
    try:
        if threads is not None:
            conn.execute(f"SET threads = {threads}")
        conn.execute(SCHEMA_PATH.read_text())

        t0 = time.monotonic()
        n_pockets, n_structures = _load_pocket_metadata(conn, screen_results, limit, scan_workers)
        stage(f"scanned {n_pockets} pockets, {n_structures} structures", t0)

        t0 = time.monotonic()
        _load_hits(conn, screen_results)
        stage("loaded leader.csv hits", t0)

        t0 = time.monotonic()
        _fill_rollups(conn)
        stage("computed rollups", t0)
    finally:
        conn.close()
    return database


def connect(database: str | os.PathLike[str], *, read_only: bool = True):
    """Open an already-built GenomeScreen database."""
    return duckdb.connect(str(database), read_only=read_only)


def summarize(database: str | os.PathLike[str]) -> str:
    """Return a short human-readable summary of a built database."""
    conn = connect(database)
    try:
        rows = [
            (
                "targets (UniProt accessions)",
                conn.sql("SELECT count(*) FROM targets").fetchone()[0],
            ),
            ("pockets", conn.sql("SELECT count(*) FROM pockets").fetchone()[0]),
            ("pocket structures", conn.sql("SELECT count(*) FROM pocket_structures").fetchone()[0]),
            ("molecules (oid x smiles)", conn.sql("SELECT count(*) FROM molecules").fetchone()[0]),
            (
                "distinct compounds",
                conn.sql("SELECT count(DISTINCT catalog_id) FROM molecules").fetchone()[0],
            ),
            ("hits", conn.sql("SELECT count(*) FROM hits").fetchone()[0]),
        ]
        by_source = conn.sql(
            "SELECT source, count(*) FROM molecules GROUP BY 1 ORDER BY 2 DESC"
        ).fetchall()
    finally:
        conn.close()
    lines = [f"{label:<28} {value:>12,}" for label, value in rows]
    lines += [f"  {source:<26} {count:>12,}" for source, count in by_source]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("root", help="GenomeScreen release directory")
    parser.add_argument("database", help="DuckDB file to create")
    parser.add_argument("--limit", type=int, help="ingest only the first N pockets")
    parser.add_argument("--overwrite", action="store_true", help="replace an existing database")
    parser.add_argument("--threads", type=int, help="DuckDB thread count")
    parser.add_argument(
        "--scan-workers",
        type=int,
        default=DEFAULT_SCAN_WORKERS,
        help=f"threads for the directory walk (default: {DEFAULT_SCAN_WORKERS})",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")
    args = parser.parse_args(argv)

    started = time.monotonic()
    path = build_database(
        args.root,
        args.database,
        limit=args.limit,
        overwrite=args.overwrite,
        threads=args.threads,
        scan_workers=args.scan_workers,
        progress=not args.quiet,
    )
    if not args.quiet:
        size_gb = path.stat().st_size / 1e9
        print(f"\n{summarize(path)}", file=sys.stderr)
        print(
            f"\nwrote {path} ({size_gb:.2f} GB) in {time.monotonic() - started:.1f}s",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
