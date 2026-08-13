"""Assemble every Parkinson's evidence layer into one DuckDB file.

Six sources, three of them independently produced, land in one database:

    PrimeKG        pd_tree, pd_edges, pathway_enrichment
    toxin overlay  toxin_target, toxin_pathway        (scripts/build_toxin_targets.py)
    ChEMBL/AACT    toxin_target_curated, target_drug, trials
    curation       mechanism_direction                (mechanism.py)
    GenomeScreen   screen_hits                        (DrugCLIP virtual screen)
    tuning         score_weights

The two toxin curations disagree, and the disagreement is kept rather than
resolved: ``curation_conflicts`` lists targets one pipeline dropped and the other
did not. The consequential one is complex I -- the upstream file discards
``nadh dehydrogenase [EC:1.6.99.3]`` as a cross-species artifact, but the human
NDUF subunits are real and rotenone's inhibition of them is the canonical PD
mechanism.

Usage:
    pdkb-build [DB_PATH] [--data-dir DIR] [--screen-source URL|none]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import duckdb

from biopharma_hackathon.pdkb import mechanism

SCHEMA_PATH = Path(__file__).with_name("schema.sql")
DEFAULT_DB = Path("pd_kb.duckdb")
DEFAULT_DATA = Path("data")
# The published GenomeScreen tables; DuckDB reads them over hf:// with predicate
# pushdown, so only the PD targets' hits are actually fetched.
DEFAULT_SCREEN = "hf://datasets/conradry/biopharma-hackathon"

NO_DRUG = "(no approved direct-MoA drug)"

# Condition-text matching for the trial rollups. Deliberately narrow: the trials
# table spans every indication these drugs were ever tested for, and the point of
# these flags is to find the small PD/neuro slice.
PD_PATTERNS = ("%parkinson%",)
NEURO_PATTERNS = (
    "%parkinson%",
    "%alzheimer%",
    "%dementia%",
    "%neurodegener%",
    "%lewy%",
    "%multiple system atrophy%",
    "%supranuclear%",
    "%huntington%",
    "%amyotrophic%",
    "%tremor%",
    "%dyskinesia%",
    "%restless legs%",
)


def _ilike_any(column: str, patterns: tuple[str, ...]) -> str:
    return " OR ".join(f"{column} ILIKE '{p}'" for p in patterns)


def _require(paths: dict[str, Path]) -> None:
    missing = {name: p for name, p in paths.items() if not p.exists()}
    if missing:
        listing = "\n".join(f"  {name}: {p}" for name, p in missing.items())
        raise SystemExit(f"missing input files:\n{listing}")


def build_database(
    database: str | os.PathLike[str] = DEFAULT_DB,
    *,
    data_dir: str | os.PathLike[str] = DEFAULT_DATA,
    screen_source: str | None = DEFAULT_SCREEN,
    overwrite: bool = True,
) -> Path:
    """Build the knowledge base and return its path."""
    db_path = Path(database)
    data = Path(data_dir)
    if overwrite and db_path.exists():
        db_path.unlink()

    sources = {
        "pd_tree": data / "pd_tree.parquet",
        "pd_subgraph": data / "pd_subgraph.parquet",
        "enrichment": data / "pd_subgraph_enrichment.parquet",
        "toxin_target": data / "pd_toxin_target.parquet",
        "toxin_pathway": data / "pd_toxin_pathway.parquet",
        "toxin_map": data / "PD_toxin_target_Toxin_Target_map.csv",
        "target_drug": data / "PD_toxin_target_Target_Drug_map.csv",
        "drug_summary": data / "PD_toxin_target_Drug_Summary.csv",
        "trials": data / "PD_toxin_target_Trials.csv",
    }
    _require(sources)

    conn = duckdb.connect(str(db_path))
    log = lambda msg: print(msg, file=sys.stderr)  # noqa: E731

    # --- PrimeKG layer -----------------------------------------------------
    conn.execute(f"""
        CREATE TABLE pd_tree AS SELECT * FROM read_parquet('{sources["pd_tree"]}');
        CREATE TABLE pd_edges AS SELECT * FROM read_parquet('{sources["pd_subgraph"]}');
        CREATE TABLE pathway_enrichment AS
            SELECT * FROM read_parquet('{sources["enrichment"]}');
    """)

    # --- toxin overlay -----------------------------------------------------
    conn.execute(f"""
        CREATE TABLE toxin_target AS SELECT * FROM read_parquet('{sources["toxin_target"]}');
        CREATE TABLE toxin_pathway AS SELECT * FROM read_parquet('{sources["toxin_pathway"]}');
    """)

    # --- the upstream ChEMBL/ClinicalTrials.gov layer -----------------------
    conn.execute(f"""
        CREATE TABLE toxin_target_curated AS
        SELECT src.Toxin                                   AS toxin,
               src.Has_human_data = 'Yes'                  AS has_human_data,
               src.N_druggable_targets                     AS n_druggable_targets,
               src.Druggable_human_targets                 AS druggable_targets,
               src.Excluded_ortholog_artifacts             AS excluded_artifacts
        FROM read_csv('{sources["toxin_map"]}', header = true, null_padding = true) src;
    """)

    # has_approved_drug separates real drug rows from the placeholder the source
    # uses to record "this target is undrugged", which is a finding, not a drug.
    conn.execute(f"""
        CREATE TABLE target_drug AS
        SELECT src.Target_gene                             AS gene,
               src.UniProt                                 AS uniprot,
               src.Protein                                 AS protein,
               src.Toxins_implicating_target               AS toxins,
               src.Drug                                    AS drug,
               try_cast(src.Approval_year AS INTEGER)      AS approval_year,
               src.Action_type                             AS action_type,
               src.Drug <> '{NO_DRUG}'                     AS has_approved_drug
        FROM read_csv('{sources["target_drug"]}',
                      header = true, nullstr = 'NULL', null_padding = true) src;
    """)

    conn.execute(f"""
        CREATE TABLE drug_summary AS
        SELECT src.Drug AS drug, src.Targets AS targets, src.Toxins AS toxins,
               try_cast(src.Approval_year AS INTEGER) AS approval_year,
               src.Action_type AS action_type,
               src.N_trials AS n_trials, src.N_biomarker_trials AS n_biomarker_trials,
               src.N_PDneuro_trials AS n_pdneuro_trials_upstream
        FROM read_csv('{sources["drug_summary"]}',
                      header = true, nullstr = 'NULL', null_padding = true) src;
    """)

    conn.execute(f"""
        CREATE TABLE trials AS
        SELECT src.Drug                                   AS drug,
               src.Targets                                AS targets,
               src.Toxins                                 AS toxins,
               src.Action                                 AS action_type,
               src.NCT_ID                                 AS nct_id,
               src.URL                                    AS url,
               src.Phase                                  AS phase,
               src.Status                                 AS status,
               try_cast(src.Start_date AS DATE)           AS start_date,
               try_cast(src.Enrollment AS INTEGER)        AS enrollment,
               src.Conditions                             AS conditions,
               try_cast(src.N_primary AS INTEGER)         AS n_primary,
               try_cast(src.N_secondary AS INTEGER)       AS n_secondary,
               src.Primary_outcome_measures               AS primary_outcomes,
               src.Secondary_outcome_measures             AS secondary_outcomes,
               coalesce(src.Biomarker_measures = 'Yes', FALSE)          AS biomarker_measures,
               coalesce({_ilike_any("src.Conditions", PD_PATTERNS)}, FALSE)    AS is_parkinsons,
               coalesce({_ilike_any("src.Conditions", NEURO_PATTERNS)}, FALSE) AS is_neuro
        FROM read_csv('{sources["trials"]}',
                      header = true, sample_size = -1, null_padding = true) src;
    """)

    # --- curated mechanism direction ---------------------------------------
    conn.execute("""
        CREATE TABLE mechanism_direction (
            gene               VARCHAR PRIMARY KEY,
            protective_actions VARCHAR,
            risk_actions       VARCHAR,
            confidence         VARCHAR NOT NULL,
            rationale          VARCHAR NOT NULL
        );
    """)
    conn.executemany("INSERT INTO mechanism_direction VALUES (?, ?, ?, ?, ?)", mechanism.rows())

    # --- GenomeScreen bridge ------------------------------------------------
    conn.execute("""
        CREATE TABLE screen_hits (
            gene VARCHAR, uniprot VARCHAR, pocket_key VARCHAR,
            catalog_id VARCHAR, source VARCHAR, smiles VARCHAR, score DOUBLE
        );
    """)
    if screen_source:
        _load_screen_hits(conn, screen_source, log)
    else:
        log("screen_hits: skipped (--screen-source none)")

    # --- scoring weights ----------------------------------------------------
    conn.execute("""
        CREATE TABLE score_weights (
            component VARCHAR PRIMARY KEY,
            weight    DOUBLE  NOT NULL,
            note      VARCHAR NOT NULL
        );
        INSERT INTO score_weights VALUES
          ('direction', 0.30, 'Does the drug oppose or mimic the toxin insult? The gate.'),
          ('pathway',   0.18, 'Enrichment of the drug target''s pathways for PD proteins.'),
          ('toxin',     0.18, 'Toxicological support for the target, specificity-weighted.'),
          ('cns',       0.17, 'Can it reach the brain? Neuro-trial precedent, with a heavy '
                              'discount for antibodies and other large molecules.'),
          ('clinical',  0.10, 'Trial volume: a proxy for accumulated human safety data. '
                              'Kept small -- it mostly measures commercial success.'),
          ('biomarker', 0.07, 'Share of trials with a biomarker readout -- measurability.');
    """)

    # --- where the two toxin curations disagree -----------------------------
    conn.execute("""
        CREATE TABLE curation_conflicts AS
        SELECT t.gene_symbol                       AS gene,
               string_agg(DISTINCT t.toxin, '; ')  AS toxins,
               max(t.mapping_confidence)           AS mapping_confidence,
               bool_or(t.in_pd_tree)               AS in_pd_tree,
               'present in the toxin overlay, absent from the upstream curated target list'
                                                   AS note
        FROM toxin_target t
        WHERE t.gene_symbol NOT IN (SELECT gene FROM target_drug)
        GROUP BY 1;
    """)

    conn.execute(SCHEMA_PATH.read_text())

    for table in (
        "pd_tree",
        "pd_edges",
        "pathway_enrichment",
        "toxin_target",
        "toxin_pathway",
        "toxin_target_curated",
        "target_drug",
        "drug_summary",
        "trials",
        "mechanism_direction",
        "screen_hits",
        "curation_conflicts",
    ):
        n = conn.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
        log(f"  {table:24} {n:>8,}")

    conn.close()
    log(f"-> {db_path}")
    return db_path


def _load_screen_hits(conn: duckdb.DuckDBPyConnection, base: str, log) -> None:
    """Pull DrugCLIP hits for the PD targets only.

    GenomeScreen is keyed by UniProt and PrimeKG by gene symbol; target_drug
    carries both, so it is the bridge. Only 17 of the 30 targets were screened at
    all -- the coverage gap is reported rather than hidden, because the targets
    that most need novel chemistry (SNCA, SOD1, SOD2) are among the missing.
    """
    try:
        conn.execute("INSTALL httpfs; LOAD httpfs;")
        conn.execute(f"""
            INSERT INTO screen_hits
            SELECT td.gene, td.uniprot, h.pocket_key, m.catalog_id, m.source, m.smiles, h.score
            FROM '{base}/hits.parquet' h
            JOIN '{base}/pockets.parquet'   p USING (pocket_key)
            JOIN '{base}/molecules.parquet' m USING (mol_id)
            JOIN (SELECT DISTINCT gene, uniprot FROM target_drug) td
                 ON td.uniprot = p.uniprot_acc;
        """)
        covered, total = conn.execute("""
            SELECT (SELECT count(DISTINCT gene) FROM screen_hits),
                   (SELECT count(DISTINCT gene) FROM target_drug)
        """).fetchone()
        log(f"screen_hits: {covered}/{total} PD targets screened by GenomeScreen")
    except Exception as exc:  # network, auth, or a moved dataset
        log(f"screen_hits: skipped ({type(exc).__name__}: {exc})")


def connect(database: str | os.PathLike[str] = DEFAULT_DB, *, read_only: bool = True):
    """Open the knowledge base."""
    return duckdb.connect(str(database), read_only=read_only)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the PD repurposing knowledge base.")
    parser.add_argument("database", nargs="?", default=str(DEFAULT_DB))
    parser.add_argument("--data-dir", default=str(DEFAULT_DATA))
    parser.add_argument(
        "--screen-source",
        default=DEFAULT_SCREEN,
        help="GenomeScreen parquet base (hf:// or local dir); 'none' to skip.",
    )
    args = parser.parse_args(argv)
    source = None if args.screen_source.lower() == "none" else args.screen_source
    build_database(args.database, data_dir=args.data_dir, screen_source=source)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
