"""Extract a Parkinson's disease subgraph from the PrimeKG edge list.

PrimeKG has no ``disease_pathway`` relation, so the disease is bridged to
pathways through its associated proteins:

    Parkinson disease
      --disease_protein-->  seed proteins
      --pathway_protein-->  pathways
      --pathway_pathway-->  (edges among those pathways)
      --pathway_protein-->  all proteins in those pathways
      --protein_protein-->  (edges among those proteins)

Bridging through proteins pulls in generic pathways alongside PD-relevant ones,
so each pathway is scored by hypergeometric enrichment of the seed proteins
against the background of all pathway-annotated proteins in the graph.

Two streaming passes over the source file: the first resolves the node sets and
scores the pathways, the second writes out every original row whose endpoints
fall inside them.

Writes two parquet files: the edge list, and a per-pathway enrichment table
(``<out>_enrichment.parquet``) ranked by q-value.

Usage:
    python scripts/extract_pd_subgraph.py [KG_CSV] [OUT_PARQUET]
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from math import exp, lgamma
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

DEFAULT_KG = Path("/home/ubuntu/kg_t1.csv")
DEFAULT_OUT = Path("data/pd_subgraph.parquet")

# x_index of the merged MONDO/Orphanet "Parkinson disease" node. The other node
# of the same name (91594, MONDO:0005180) carries only disease_disease edges.
PD_INDEX = "80275"

HEADER = [
    "relation",
    "display_relation",
    "x_index",
    "x_id",
    "x_type",
    "x_name",
    "x_source",
    "y_index",
    "y_id",
    "y_type",
    "y_name",
    "y_source",
]


def endpoints(row: list[str]) -> tuple[str, str]:
    return row[2], row[7]


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(N population, K successes, n draws).

    Summed in log space so the factorials stay finite at PrimeKG's scale.
    """
    if k <= 0:
        return 1.0
    denom = _log_comb(N, n)
    total = 0.0
    for i in range(k, min(n, K) + 1):
        term = _log_comb(K, i) + _log_comb(N - K, n - i) - denom
        if term > float("-inf"):
            total += exp(term)
    return min(total, 1.0)


def benjamini_hochberg(pvalues: dict[str, float]) -> dict[str, float]:
    """BH step-up FDR correction, returning monotonic q-values."""
    ordered = sorted(pvalues.items(), key=lambda kv: kv[1])
    m = len(ordered)
    qvalues: dict[str, float] = {}
    running = 1.0
    for rank in range(m, 0, -1):
        key, p = ordered[rank - 1]
        running = min(running, p * m / rank, 1.0)
        qvalues[key] = running
    return qvalues


@dataclass
class Scan:
    """Everything pass 1 learns about the Parkinson's neighbourhood."""

    seed_proteins: set[str]
    pathways: set[str]
    proteins: set[str]
    stats: dict[str, dict[str, object]]
    pathway_members: dict[str, set[str]]
    protein_names: dict[str, str]
    # protein index -> {(drug index, display_relation): (drug_id, drug_name)}
    drugs_by_protein: dict[str, dict[tuple[str, str], tuple[str, str]]]


def resolve_nodes(kg_path: Path, *, with_drugs: bool = False) -> Scan:
    """Pass 1: find seed proteins, their pathways, those pathways' proteins, and
    the hypergeometric enrichment of each pathway for the seed proteins.

    Set ``with_drugs`` to also collect drug_protein edges for those proteins.
    PrimeKG has no drug_pathway relation, so proteins are the only drug bridge.
    """
    disease_protein: list[tuple[str, str]] = []
    pathway_members: dict[str, set[str]] = {}
    pathway_names: dict[str, str] = {}
    protein_names: dict[str, str] = {}
    universe: set[str] = set()
    drug_edges: list[tuple[str, str, str, str]] = []

    with kg_path.open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            rel = row[0]
            if rel == "disease_protein":
                disease_protein.append(endpoints(row))
            elif rel == "pathway_protein":
                x, y = endpoints(row)
                if row[4] == "pathway":
                    pathway, protein, name, pname = x, y, row[5], row[10]
                else:
                    pathway, protein, name, pname = y, x, row[10], row[5]
                pathway_members.setdefault(pathway, set()).add(protein)
                pathway_names[pathway] = name
                protein_names[protein] = pname
                universe.add(protein)
            elif with_drugs and rel == "drug_protein":
                x, y = endpoints(row)
                if row[4] == "drug":
                    drug, protein, drug_id, drug_name = x, y, row[3], row[5]
                else:
                    drug, protein, drug_id, drug_name = y, x, row[8], row[10]
                # Stored in both orientations; the dict collapses the duplicate.
                drug_edges.append((protein, drug, drug_id, drug_name, row[1]))

    seed_proteins = {y if x == PD_INDEX else x for x, y in disease_protein if PD_INDEX in (x, y)}

    pathways = {pathway for pathway, members in pathway_members.items() if members & seed_proteins}
    proteins = set().union(*(pathway_members[p] for p in pathways))

    # Only seed proteins that are pathway-annotated can be drawn, so they define
    # the sample size; the rest carry no information about pathway membership.
    annotated_seed = seed_proteins & universe
    N, n = len(universe), len(annotated_seed)

    stats: dict[str, dict[str, object]] = {}
    for pathway in pathways:
        members = pathway_members[pathway]
        K = len(members)
        k = len(members & annotated_seed)
        expected = n * K / N
        stats[pathway] = {
            "pathway_name": pathway_names[pathway],
            "pathway_size": K,
            "seed_overlap": k,
            "expected_overlap": round(expected, 4),
            "fold_enrichment": round(k / expected, 4) if expected else "",
            "pvalue": hypergeom_sf(k, N, K, n),
        }

    qvalues = benjamini_hochberg({p: s["pvalue"] for p, s in stats.items()})
    for pathway, s in stats.items():
        s["qvalue"] = qvalues[pathway]

    # Keyed by (drug, display_relation), not drug alone: 185 drug-protein pairs
    # carry two relation types (usually enzyme + target) and both are real edges.
    drugs_by_protein: dict[str, dict[tuple[str, str], tuple[str, str]]] = {}
    for protein, drug, drug_id, drug_name, display in drug_edges:
        if protein in proteins:
            drugs_by_protein.setdefault(protein, {})[(drug, display)] = (drug_id, drug_name)

    print(
        f"enrichment background: {N} pathway-annotated proteins, "
        f"{n}/{len(seed_proteins)} seed proteins annotated",
        file=sys.stderr,
    )
    return Scan(
        seed_proteins=seed_proteins,
        pathways=pathways,
        proteins=proteins,
        stats=stats,
        pathway_members=pathway_members,
        protein_names=protein_names,
        drugs_by_protein=drugs_by_protein,
    )


def keep(row: list[str], pathways: set[str], proteins: set[str]) -> str | None:
    """Return the subgraph layer this row belongs to, or None to drop it."""
    rel = row[0]
    x, y = endpoints(row)

    if rel == "disease_protein" and PD_INDEX in (x, y):
        return "disease_protein"
    if rel == "pathway_pathway" and x in pathways and y in pathways:
        return "pathway_pathway"
    if rel == "pathway_protein":
        pathway, protein = (x, y) if row[4] == "pathway" else (y, x)
        if pathway in pathways and protein in proteins:
            return "pathway_protein"
    if rel == "protein_protein" and x in proteins and y in proteins:
        return "protein_protein"
    return None


def main() -> None:
    kg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_KG
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUT
    enrichment_path = out_path.with_name(f"{out_path.stem}_enrichment.parquet")

    scan = resolve_nodes(kg_path)
    pathways, proteins, stats = scan.pathways, scan.proteins, scan.stats
    print(
        f"seed proteins: {len(scan.seed_proteins)}  "
        f"pathways: {len(pathways)}  "
        f"pathway proteins: {len(proteins)}",
        file=sys.stderr,
    )

    columns: dict[str, list] = {name: [] for name in HEADER}
    columns["layer"] = []
    columns["x_pathway_qvalue"] = []
    columns["y_pathway_qvalue"] = []
    counts: dict[str, int] = {}

    with kg_path.open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            layer = keep(row, pathways, proteins)
            if layer is None:
                continue
            for name, value in zip(HEADER, row, strict=True):
                columns[name].append(int(value) if name in ("x_index", "y_index") else value)
            columns["layer"].append(layer)
            x, y = endpoints(row)
            columns["x_pathway_qvalue"].append(stats[x]["qvalue"] if row[4] == "pathway" else None)
            columns["y_pathway_qvalue"].append(stats[y]["qvalue"] if row[9] == "pathway" else None)
            counts[layer] = counts.get(layer, 0) + 1

    schema = pa.schema(
        [
            ("relation", pa.string()),
            ("display_relation", pa.string()),
            ("x_index", pa.int32()),
            ("x_id", pa.string()),
            ("x_type", pa.string()),
            ("x_name", pa.string()),
            ("x_source", pa.string()),
            ("y_index", pa.int32()),
            ("y_id", pa.string()),
            ("y_type", pa.string()),
            ("y_name", pa.string()),
            ("y_source", pa.string()),
            ("layer", pa.string()),
            ("x_pathway_qvalue", pa.float64()),
            ("y_pathway_qvalue", pa.float64()),
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    edges = pa.table(
        {name: pa.array(values, type=schema.field(name).type) for name, values in columns.items()},
        schema=schema,
    )
    pq.write_table(edges, out_path, compression="zstd")

    ranked = sorted(stats.items(), key=lambda kv: (kv[1]["qvalue"], kv[1]["pvalue"]))
    enrichment = pa.table(
        {
            "pathway_index": pa.array([int(p) for p, _ in ranked], pa.int32()),
            "pathway_name": pa.array([s["pathway_name"] for _, s in ranked], pa.string()),
            "pathway_size": pa.array([s["pathway_size"] for _, s in ranked], pa.int32()),
            "seed_overlap": pa.array([s["seed_overlap"] for _, s in ranked], pa.int32()),
            "expected_overlap": pa.array([s["expected_overlap"] for _, s in ranked], pa.float64()),
            "fold_enrichment": pa.array(
                [s["fold_enrichment"] or None for _, s in ranked], pa.float64()
            ),
            "pvalue": pa.array([s["pvalue"] for _, s in ranked], pa.float64()),
            "qvalue": pa.array([s["qvalue"] for _, s in ranked], pa.float64()),
        }
    )
    pq.write_table(enrichment, enrichment_path, compression="zstd")

    for layer in ("disease_protein", "pathway_pathway", "pathway_protein", "protein_protein"):
        print(f"{layer}: {counts.get(layer, 0)}", file=sys.stderr)
    print(f"total: {sum(counts.values())} -> {out_path}", file=sys.stderr)
    significant = sum(1 for _, s in ranked if s["qvalue"] < 0.05)
    print(f"pathways at q<0.05: {significant}/{len(ranked)} -> {enrichment_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
