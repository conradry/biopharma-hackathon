"""Tests for the PD knowledge base, run against a synthetic mini-corpus.

The fixture is small but carries one of each case the scoring has to get right:
a drug that opposes the toxin insult (MAO-B inhibitor), one that mimics it (AKT
inhibitor), one whose direction is genuinely two-sided (cholinesterase), an
antibody that should be discounted for a CNS indication, and a target with no
approved drug at all.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biopharma_hackathon.pdkb import build_database, connect, evidence_graph
from biopharma_hackathon.pdkb.mechanism import classify

PATHWAYS = {
    # index: (name, size, qvalue)
    1: ("Biogenic amines are oxidatively deaminated by MAOA and MAOB", 2, 0.0064),
    2: ("Interleukin-4 and Interleukin-13 signaling", 108, 0.0064),
    3: ("Neurotransmitter clearance", 4, 0.117),
}
# gene -> (pathway index, is_pd_seed)
GENES = {
    "MAOB": (1, True),
    "AKT1": (2, False),
    "ACHE": (3, False),
    "TNF": (2, True),
    "SNCA": (2, True),
}
# drug -> (gene, action, approval year)
DRUGS = {
    "Rasagiline": ("MAOB", "INHIBITOR", 2006),
    "Capivasertib": ("AKT1", "INHIBITOR", 2023),
    "Donepezil": ("ACHE", "INHIBITOR", 1996),
    "Adalimumab": ("TNF", "INHIBITOR", 2002),
}


def _write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


@pytest.fixture(scope="module")
def kb(tmp_path_factory) -> Path:
    data = tmp_path_factory.mktemp("data")

    tree = []
    for gene, (pathway, is_seed) in GENES.items():
        name, size, qvalue = PATHWAYS[pathway]
        base = {
            "disease_index": 80275,
            "disease_name": "Parkinson disease",
            "pathway_index": pathway,
            "pathway_name": name,
            "pathway_size": size,
            "pathway_seed_overlap": 2,
            "pathway_fold_enrichment": 12.0,
            "pathway_qvalue": qvalue,
            "protein_index": abs(hash(gene)) % 10000,
            "protein_name": gene,
            "protein_is_pd_seed": is_seed,
        }
        drug = next((d for d, (g, *_) in DRUGS.items() if g == gene), None)
        tree.append(
            {
                **base,
                "drug_index": abs(hash(drug)) % 10000 if drug else None,
                "drug_id": "DB00001" if drug else None,
                "drug_name": drug,
                "drug_relation": "target" if drug else None,
                "depth": 4 if drug else 3,
            }
        )
    pq.write_table(pa.Table.from_pylist(tree), data / "pd_tree.parquet")
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "relation": "pathway_protein",
                    "x_index": 1,
                    "y_index": 2,
                    "layer": "pathway_protein",
                }
            ]
        ),
        data / "pd_subgraph.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "pathway_index": i,
                    "pathway_name": n,
                    "pathway_size": s,
                    "seed_overlap": 2,
                    "expected_overlap": 0.1,
                    "fold_enrichment": 12.0,
                    "pvalue": q / 2,
                    "qvalue": q,
                }
                for i, (n, s, q) in PATHWAYS.items()
            ]
        ),
        data / "pd_subgraph_enrichment.parquet",
    )

    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "toxin": "mptp",
                    "gene_symbol": "MAOB",
                    "specificity_weight": 0.09,
                    "in_pd_tree": True,
                    "mapping_confidence": "family",
                },
                {
                    "toxin": "rotenone",
                    "gene_symbol": "AKT1",
                    "specificity_weight": 0.20,
                    "in_pd_tree": True,
                    "mapping_confidence": "exact",
                },
                {
                    "toxin": "paraquat",
                    "gene_symbol": "TNF",
                    "specificity_weight": 0.35,
                    "in_pd_tree": True,
                    "mapping_confidence": "exact",
                },
                {
                    "toxin": "rotenone",
                    "gene_symbol": "NDUFS1",
                    "specificity_weight": 0.11,
                    "in_pd_tree": False,
                    "mapping_confidence": "complex",
                },
            ]
        ),
        data / "pd_toxin_target.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "pathway_index": 2,
                    "n_toxins": 2,
                    "n_toxin_targets": 2,
                    "toxin_support": 0.55,
                    "toxins": "paraquat; rotenone",
                    "toxin_targets": "AKT1; TNF",
                }
            ]
        ),
        data / "pd_toxin_pathway.parquet",
    )

    _write_csv(
        data / "PD_toxin_target_Toxin_Target_map.csv",
        [
            "Toxin",
            "Has_human_data",
            "N_druggable_targets",
            "Druggable_human_targets",
            "Excluded_ortholog_artifacts",
        ],
        [["mptp", "Yes", 1, "MAOB"], ["rotenone", "Yes", 1, "AKT1", "nadh dehydrogenase"]],
    )
    _write_csv(
        data / "PD_toxin_target_Target_Drug_map.csv",
        [
            "Target_gene",
            "UniProt",
            "Protein",
            "Toxins_implicating_target",
            "Drug",
            "Approval_year",
            "Action_type",
        ],
        [
            [gene, "P00000", f"{gene} protein", "mptp", drug, year, action]
            for drug, (gene, action, year) in DRUGS.items()
        ]
        # a target nothing approved engages: must not become a candidate row
        + [
            [
                "SNCA",
                "P37840",
                "Alpha-synuclein",
                "rotenone",
                "(no approved direct-MoA drug)",
                "NULL",
                "NULL",
            ]
        ],
    )
    _write_csv(
        data / "PD_toxin_target_Drug_Summary.csv",
        [
            "Drug",
            "Targets",
            "Toxins",
            "Approval_year",
            "Action_type",
            "N_trials",
            "N_biomarker_trials",
            "N_PDneuro_trials",
        ],
        [[d, g, "mptp", y, a, 2, 1, 1] for d, (g, a, y) in DRUGS.items()],
    )
    trials = []
    for drug in DRUGS:
        condition = "Parkinson's Disease" if drug in ("Rasagiline", "Donepezil") else "Cancer"
        for i in range(2):
            trials.append(
                [
                    "mptp",
                    DRUGS[drug][0],
                    drug,
                    DRUGS[drug][2],
                    DRUGS[drug][1],
                    f"NCT{abs(hash(drug + str(i))) % 10**8:08d}",
                    "https://example.org",
                    "Phase 3",
                    "Completed",
                    "2010-01-01",
                    100,
                    condition,
                    1,
                    1,
                    "primary measure",
                    "secondary measure",
                    "Yes" if i else "No",
                ]
            )
    _write_csv(
        data / "PD_toxin_target_Trials.csv",
        [
            "Toxins",
            "Targets",
            "Drug",
            "Approval_yr",
            "Action",
            "NCT_ID",
            "URL",
            "Phase",
            "Status",
            "Start_date",
            "Enrollment",
            "Conditions",
            "N_primary",
            "N_secondary",
            "Primary_outcome_measures",
            "Secondary_outcome_measures",
            "Biomarker_measures",
        ],
        trials,
    )

    db = data / "kb.duckdb"
    build_database(db, data_dir=data, screen_source=None)
    return db


@pytest.fixture(scope="module")
def conn(kb: Path):
    with connect(kb) as connection:
        yield connection


def _candidate(conn, drug: str) -> dict:
    row = conn.execute("SELECT * FROM drug_candidates WHERE drug = ?", [drug]).fetchone()
    names = [d[0] for d in conn.description]
    return dict(zip(names, row, strict=True))


def test_direction_gates_the_verdict(conn) -> None:
    """The whole point: same evidence, opposite direction, opposite conclusion."""
    assert _candidate(conn, "Rasagiline")["verdict"] == "candidate"
    assert _candidate(conn, "Capivasertib")["verdict"] == "mechanistic_risk"


def test_ambiguous_direction_is_not_scored_as_support(conn) -> None:
    donepezil = _candidate(conn, "Donepezil")
    assert donepezil["n_protective"] == 0
    assert donepezil["n_ambiguous"] == 1
    # 0.5 is the neutral midpoint: no directional information either way
    assert donepezil["direction_component"] == pytest.approx(0.5)


def test_risk_drug_scores_below_protective_drug(conn) -> None:
    assert (
        _candidate(conn, "Capivasertib")["priority_score"]
        < _candidate(conn, "Rasagiline")["priority_score"]
    )


def test_antibody_is_discounted_for_a_cns_indication(conn) -> None:
    adalimumab = _candidate(conn, "Adalimumab")
    assert adalimumab["modality"] == "biologic"
    # no PD and no neuro trials in the fixture, so it takes the biologic floor
    assert adalimumab["cns_component"] == pytest.approx(0.15)


def test_undrugged_target_is_not_a_candidate(conn) -> None:
    drugs = {r[0] for r in conn.execute("SELECT drug FROM drug_candidates").fetchall()}
    assert "(no approved direct-MoA drug)" not in drugs
    undrugged = {r[0] for r in conn.execute("SELECT gene FROM undrugged_targets").fetchall()}
    assert undrugged == {"SNCA"}


def test_curation_conflict_records_complex_i(conn) -> None:
    """The upstream list drops complex I; the overlay keeps it. Keep the disagreement."""
    conflicts = {r[0] for r in conn.execute("SELECT gene FROM curation_conflicts").fetchall()}
    assert "NDUFS1" in conflicts


def test_weights_are_tunable_without_a_rebuild(kb: Path, tmp_path: Path) -> None:
    # A copy, so the module-scoped read-only connection stays open alongside it.
    writable = tmp_path / "tunable.duckdb"
    writable.write_bytes(kb.read_bytes())
    with connect(writable, read_only=False) as conn:
        before = _candidate(conn, "Rasagiline")["priority_score"]
        conn.execute("UPDATE score_weights SET weight = 0 WHERE component = 'direction'")
        after = _candidate(conn, "Rasagiline")["priority_score"]
    assert after < before


def test_evidence_graph_renders_in_every_format(conn) -> None:
    graph = evidence_graph(conn, "Rasagiline")
    kinds = {n.kind for n in graph.nodes}
    assert {"drug", "target", "toxin", "pathway", "disease"} <= kinds
    assert graph.summary["verdict"] == "candidate"

    mermaid, dot, blob = graph.to_mermaid(), graph.to_dot(), graph.to_json()
    assert mermaid.startswith("flowchart LR")
    assert "class n1 protective;" in mermaid
    assert dot.startswith('digraph "Rasagiline"') and dot.rstrip().endswith("}")
    assert '"MAOB"' in blob


def test_evidence_graph_marks_a_risk_edge_red(conn) -> None:
    graph = evidence_graph(conn, "Capivasertib")
    acts_on = [e for e in graph.edges if e.kind == "acts_on"]
    assert [e.attrs["direction"] for e in acts_on] == ["risk"]
    assert "risk;" in graph.to_mermaid()
    assert "#cf222e" in graph.to_dot()


def test_graph_lookup_is_case_insensitive_and_errors_cleanly(conn) -> None:
    assert evidence_graph(conn, "rasagiline").drug == "Rasagiline"
    with pytest.raises(LookupError):
        evidence_graph(conn, "not-a-drug")


@pytest.mark.parametrize(
    ("gene", "action", "expected"),
    [
        ("MAOB", "INHIBITOR", "protective"),
        ("MAOB", "AGONIST", "risk"),
        ("AKT1", "INHIBITOR", "risk"),
        ("AKT1", "ACTIVATOR", "protective"),
        ("SLC6A3", "RELEASING AGENT", "risk"),
        ("SLC6A3", "INHIBITOR", "protective"),
        ("TH", "INHIBITOR", "risk"),
        ("ACHE", "INHIBITOR", "ambiguous"),
        ("ALB", "INHIBITOR", "unknown"),
        ("MAOB", None, "unknown"),
        ("NOT_A_GENE", "INHIBITOR", "unknown"),
    ],
)
def test_classify(gene: str, action: str | None, expected: str) -> None:
    assert classify(gene, action) == expected
