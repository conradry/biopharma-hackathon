"""Tests for the GenomeScreen ingest, run against a synthetic mini-release.

The fixture reproduces the quirks that matter: template vs. detected pockets, a missing
grid file, an oid appearing under two protomer SMILES, and both vendor id namespaces.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from biopharma_hackathon.genomescreen import (
    build_database,
    connect,
    parse_grid_file,
    parse_pocket_dirname,
    scan_pocket_dir,
)
from biopharma_hackathon.genomescreen.ingest import resolve_screen_results

TEMPLATE_POCKET = "AF-P12345-F1-model_v4_0_0"
DETECTED_POCKET = "AF-P12345-F1-model_v4_1_pocket3"
OTHER_POCKET = "AF-Q9Y6K9-F1-model_v4_0_pocket1"


def _write_structure(directory: Path, stem: str, center: tuple[float, float, float] | None) -> None:
    (directory / f"{stem}.pdbgz").write_bytes(gzip.compress(b"ATOM  \nEND\n"))
    if center is None:
        return
    x, y, z = center
    (directory / f"{stem}_grid.in").write_text(
        f"GRIDFILE {stem}_grid.zip\nGRID_CENTER {x}, {y}, {z}\nRECEP_FILE {stem}.maegz\n"
    )


@pytest.fixture(scope="module")
def release(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("GenomeScreenMini")
    screen = root / "screen_results"

    template = screen / TEMPLATE_POCKET
    template.mkdir(parents=True)
    _write_structure(template, f"{TEMPLATE_POCKET}_0_5vm0_complex_refined", (1.0, 2.0, -3.5))
    _write_structure(template, f"{TEMPLATE_POCKET}_1_3qxt_complex_refined", (1.5, 2.5, -3.0))
    (template / "leader.csv").write_text(
        "smiles,Name,oid,score\n"
        "OCc1coc2cc(I)ccc12,22073,ZINC000066055208,7.83\n"
        "FC1=CC=C(C=C1Cl)C1=NNC=C1,78967,Z1333761449_1_T2,7.56\n"
        "O=c1[nH]c(=S)[nH]nc1-c1ccc(O)cc1,36,ZINC000004243184,7.10\n"
    )

    detected = screen / DETECTED_POCKET
    detected.mkdir()
    # This structure ships without a grid file, as three do in the public release.
    _write_structure(detected, f"{DETECTED_POCKET}_0_complex_refined", None)
    _write_structure(detected, f"{DETECTED_POCKET}_1_complex_refined", (0.0, 0.0, 0.0))
    (detected / "leader.csv").write_text(
        "smiles,Name,oid,score\n"
        # Same oid as above, different protomer -> a second molecules row, one catalog_id.
        "O=c1[n-]c(=S)[nH]nc1-c1ccc(O)cc1,259,ZINC000004243184,8.40\n"
        "CCO,1,PV-001914042032_1_T1,6.20\n"
    )

    other = screen / OTHER_POCKET
    other.mkdir()
    _write_structure(other, f"{OTHER_POCKET}_0_complex_refined", (9.0, 9.0, 9.0))
    (other / "leader.csv").write_text("smiles,Name,oid,score\nCCN,7,Z999_1,5.00\n")

    return root


@pytest.fixture(scope="module")
def db(release: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("db") / "genomescreen.duckdb"
    return build_database(release, path)


def test_parse_template_pocket_dirname():
    pocket = parse_pocket_dirname(TEMPLATE_POCKET)
    assert pocket.uniprot_acc == "P12345"
    assert pocket.af_model_version == 4
    assert pocket.fragment_idx == 0
    assert pocket.pocket_kind == "template"
    assert pocket.pocket_idx == 0


def test_parse_detected_pocket_dirname():
    pocket = parse_pocket_dirname(DETECTED_POCKET)
    assert (pocket.fragment_idx, pocket.pocket_kind, pocket.pocket_idx) == (1, "detected", 3)


def test_parse_pocket_dirname_rejects_junk():
    with pytest.raises(ValueError):
        parse_pocket_dirname("not-a-pocket")


def test_parse_grid_file(release: Path):
    grid = parse_grid_file(
        release
        / "screen_results"
        / TEMPLATE_POCKET
        / f"{TEMPLATE_POCKET}_0_5vm0_complex_refined_grid.in"
    )
    assert grid["grid_center"] == (1.0, 2.0, -3.5)
    assert grid["receptor_file"].endswith(".maegz")


def test_scan_pocket_dir_extracts_template_pdb_ids(release: Path):
    _, structures = scan_pocket_dir(release / "screen_results" / TEMPLATE_POCKET)
    assert [s.structure_idx for s in structures] == [0, 1]
    assert [s.template_pdb_id for s in structures] == ["5vm0", "3qxt"]


def test_scan_pocket_dir_tolerates_missing_grid(release: Path):
    _, structures = scan_pocket_dir(release / "screen_results" / DETECTED_POCKET)
    missing = next(s for s in structures if s.structure_idx == 0)
    assert missing.grid_path is None
    assert missing.grid_center_x is None
    assert missing.structure_path.endswith(".pdbgz")


def test_resolve_screen_results_accepts_either_level(release: Path):
    assert resolve_screen_results(release) == release / "screen_results"
    assert resolve_screen_results(release / "screen_results") == release / "screen_results"


def test_row_counts(db: Path):
    conn = connect(db)
    assert conn.sql("SELECT count(*) FROM targets").fetchone()[0] == 2
    assert conn.sql("SELECT count(*) FROM pockets").fetchone()[0] == 3
    assert conn.sql("SELECT count(*) FROM pocket_structures").fetchone()[0] == 5
    assert conn.sql("SELECT count(*) FROM hits").fetchone()[0] == 6
    # 6 hit rows, but ZINC000004243184 contributes two (oid, smiles) molecule rows.
    assert conn.sql("SELECT count(*) FROM molecules").fetchone()[0] == 6
    assert conn.sql("SELECT count(DISTINCT catalog_id) FROM molecules").fetchone()[0] == 5


def test_source_and_catalog_id_assignment(db: Path):
    conn = connect(db)
    rows = dict(conn.sql("SELECT DISTINCT oid, source FROM molecules").fetchall())
    assert rows["ZINC000066055208"] == "ZINC"
    assert rows["Z1333761449_1_T2"] == "Enamine REAL"
    assert rows["Z999_1"] == "Enamine REAL"
    assert rows["PV-001914042032_1_T1"] == "Enamine PV"

    catalog, protomer, tautomer = conn.sql(
        "SELECT catalog_id, protomer_idx, tautomer_idx FROM molecules WHERE oid = 'Z1333761449_1_T2'"
    ).fetchone()
    assert (catalog, protomer, tautomer) == ("Z1333761449", 1, 2)
    # ZINC ids carry no protomer/tautomer suffix.
    assert conn.sql(
        "SELECT catalog_id, protomer_idx FROM molecules WHERE oid = 'ZINC000066055208'"
    ).fetchone() == ("ZINC000066055208", None)


def test_protomer_variants_share_a_catalog_id(db: Path):
    conn = connect(db)
    smiles = conn.sql(
        "SELECT smiles FROM molecules WHERE catalog_id = 'ZINC000004243184' ORDER BY smiles"
    ).fetchall()
    assert len(smiles) == 2
    assert {s for (s,) in smiles} == {
        "O=c1[nH]c(=S)[nH]nc1-c1ccc(O)cc1",
        "O=c1[n-]c(=S)[nH]nc1-c1ccc(O)cc1",
    }


def test_rank_and_rollups(db: Path):
    conn = connect(db)
    ranked = conn.sql(
        "SELECT rank_in_pocket, score FROM hits WHERE pocket_key = ? ORDER BY rank_in_pocket",
        params=[TEMPLATE_POCKET],
    ).fetchall()
    assert ranked == [(1, 7.83), (2, 7.56), (3, 7.10)]

    assert conn.sql(
        "SELECT n_hits, n_structures, best_score FROM pockets WHERE pocket_key = ?",
        params=[TEMPLATE_POCKET],
    ).fetchone() == (3, 2, 7.83)

    assert conn.sql(
        "SELECT n_pockets, n_hits, best_score FROM targets WHERE uniprot_acc = 'P12345'"
    ).fetchone() == (2, 5, 8.40)


def test_hit_details_view_joins_through(db: Path):
    conn = connect(db)
    row = conn.sql(
        "SELECT uniprot_acc, pocket_kind, oid, source, score FROM hit_details "
        "ORDER BY score DESC LIMIT 1"
    ).fetchone()
    assert row == ("P12345", "detected", "ZINC000004243184", "ZINC", 8.40)


def test_target_molecule_best_collapses_pockets(db: Path):
    conn = connect(db)
    assert conn.sql(
        "SELECT n_hits, n_pockets, best_score, best_smiles FROM target_molecule_best "
        "WHERE uniprot_acc = 'P12345' AND catalog_id = 'ZINC000004243184'"
    ).fetchone() == (2, 2, 8.40, "O=c1[n-]c(=S)[nH]nc1-c1ccc(O)cc1")


def test_limit_keeps_hits_consistent(release: Path, tmp_path: Path):
    path = build_database(release, tmp_path / "partial.duckdb", limit=1)
    conn = connect(path)
    assert conn.sql("SELECT count(*) FROM pockets").fetchone()[0] == 1
    # Hits from the two skipped pockets must not leak in via the CSV glob.
    assert conn.sql("SELECT count(*) FROM hits").fetchone()[0] == 3
    assert conn.sql("SELECT count(*) FROM molecules").fetchone()[0] == 3


def test_overwrite_guard(release: Path, tmp_path: Path):
    path = tmp_path / "twice.duckdb"
    build_database(release, path, limit=1)
    with pytest.raises(FileExistsError):
        build_database(release, path, limit=1)
    build_database(release, path, limit=1, overwrite=True)
