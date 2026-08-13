"""Export the top-10 repurposing ranking and its evidence graphs.

Produces one ranking table and a set of graph files small enough to actually
look at. The edge budget is the design constraint: a graph of everything the
knowledge base knows is unreadable, so each export is capped (default 100 edges)
and the script reports what it dropped rather than silently truncating.

Two graphs, because they answer different questions:

    mechanism   toxin -> target <- drug, target -> pathway -> Parkinson's
                why these drugs are on the list at all
    readout     drug -> biomarker category
                how you would tell whether one is working

The ranking excludes drugs already used in PD and drugs whose direction is wrong
(see mechanism.py) -- those are validation and anti-recommendations, not
proposals. The three approved PD drugs ride along in the mechanism graph as
labelled anchors so a viewer can see the method recovering known biology.

Outputs (to --out-dir, default data/):
    pd_ranking_top10.csv        the ranking, one row per drug
    pd_top10_nodes.csv          node table (Gephi / Cytoscape)
    pd_top10_edges.csv          edge table
    pd_top10.graphml            same graph, single file
    pd_top10_mechanism.mmd      Mermaid, paste into markdown
    pd_top10_readouts.mmd       Mermaid, biomarker readouts
    pd_pathway_convergence.mmd  Mermaid, pathway-level view

Usage:
    python scripts/export_ranking_graph.py [--database pd_kb.duckdb]
                                           [--top 10] [--max-edges 100]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from xml.sax.saxutils import escape

import duckdb

DEFAULT_DB = Path("pd_kb.duckdb")
DEFAULT_OUT = Path("data")
PD_ROOT = "Parkinson disease"

RANKING_SQL = """
SELECT drug, targets, toxins, modality, verdict, approval_year,
       priority_score, direction_component, pathway_component, toxin_component,
       cns_component, clinical_component, biomarker_component,
       n_protective, n_risk, n_ambiguous,
       best_pathway_qvalue, n_enriched_pathways, toxin_support, max_toxins_per_target,
       n_trials, n_late_phase_trials, n_pd_trials, n_neuro_trials,
       n_disease_readout_trials
FROM drug_candidates
WHERE verdict = 'candidate' AND n_pd_trials = 0
ORDER BY priority_score DESC
LIMIT ?
"""

# Approved PD drugs, carried into the mechanism graph as anchors. They are not
# proposals -- the point is that the ranking found them without being told.
ANCHOR_SQL = """
SELECT drug, targets, priority_score, n_pd_trials
FROM drug_candidates
WHERE verdict = 'established_pd_therapy'
ORDER BY priority_score DESC
LIMIT 3
"""

DIRECTION_COLOUR = {
    "protective": "#1a7f37",
    "risk": "#cf222e",
    "ambiguous": "#bf8700",
    "unknown": "#6e7781",
}


class Graph:
    """Nodes and edges with a hard edge budget."""

    def __init__(self, max_edges: int) -> None:
        self.max_edges = max_edges
        self.nodes: dict[str, dict] = {}
        self.edges: list[dict] = []
        self.dropped = 0

    def node(self, node_id: str, kind: str, label: str, **attrs) -> str:
        if node_id not in self.nodes:
            self.nodes[node_id] = {"id": node_id, "kind": kind, "label": label, **attrs}
        return node_id

    def edge(self, source: str, target: str, kind: str, label: str = "", **attrs) -> bool:
        if len(self.edges) >= self.max_edges:
            self.dropped += 1
            return False
        self.edges.append(
            {"source": source, "target": target, "kind": kind, "label": label, **attrs}
        )
        return True

    def prune_orphans(self) -> None:
        """Drop nodes that lost every edge to the budget."""
        touched = {e["source"] for e in self.edges} | {e["target"] for e in self.edges}
        self.nodes = {k: v for k, v in self.nodes.items() if k in touched}


def build_mechanism_graph(conn, drugs: list[str], anchors: list[str], max_edges: int) -> Graph:
    """toxin -> target <- drug, target -> pathway -> Parkinson's."""
    graph = Graph(max_edges)
    disease = graph.node(f"disease:{PD_ROOT}", "disease", PD_ROOT)
    everyone = drugs + anchors

    rows = conn.execute(
        """
        SELECT e.drug, e.gene, e.action_type, e.direction, e.direction_confidence,
               e.toxins, e.n_toxins, e.toxin_support, e.best_pathway_qvalue,
               c.priority_score, c.verdict, c.modality, c.n_trials,
               c.n_disease_readout_trials
        FROM drug_target_evidence e
        JOIN drug_candidates c USING (drug)
        WHERE e.drug IN ?
        ORDER BY c.priority_score DESC, e.gene
        """,
        [everyone],
    ).fetchall()

    genes: set[str] = set()
    for (
        drug,
        gene,
        action,
        direction,
        confidence,
        toxins,
        n_toxins,
        support,
        best_q,
        score,
        verdict,
        modality,
        n_trials,
        n_readouts,
    ) in rows:
        is_anchor = verdict == "established_pd_therapy"
        drug_id = graph.node(
            f"drug:{drug}",
            "anchor_drug" if is_anchor else "drug",
            drug,
            priority_score=score,
            verdict=verdict,
            modality=modality,
            n_trials=n_trials,
            n_disease_readout_trials=n_readouts,
            rank=0 if is_anchor else drugs.index(drug) + 1,
        )
        gene_id = graph.node(
            f"target:{gene}",
            "target",
            gene,
            n_toxins=n_toxins or 0,
            toxin_support=round(support or 0, 4),
            best_pathway_qvalue=best_q,
        )
        genes.add(gene)
        graph.edge(
            drug_id,
            gene_id,
            "acts_on",
            (action or "unknown").lower(),
            direction=direction,
            direction_confidence=confidence or "",
        )
        for toxin in filter(None, (t.strip() for t in (toxins or "").split(";"))):
            graph.edge(
                graph.node(f"toxin:{toxin}", "toxin", toxin), gene_id, "implicates", "implicates"
            )

    # Enriched pathways only, and at most two per target -- pathway membership is
    # the bushiest part of the graph and the least informative beyond the top hits.
    for gene in sorted(genes):
        for pathway, qvalue, fold, size in conn.execute(
            """
            SELECT pathway_name, min(pathway_qvalue), max(pathway_fold_enrichment),
                   max(pathway_size)
            FROM pd_tree
            WHERE protein_name = ? AND pathway_qvalue < 0.05
            GROUP BY 1 ORDER BY 2 LIMIT 2
            """,
            [gene],
        ).fetchall():
            pathway_id = graph.node(
                f"pathway:{pathway}",
                "pathway",
                pathway,
                qvalue=qvalue,
                fold_enrichment=fold,
                pathway_size=size,
            )
            graph.edge(f"target:{gene}", pathway_id, "member_of", "in pathway")
            graph.edge(pathway_id, disease, "enriched_in", f"q={qvalue:.2g}", qvalue=qvalue)

    graph.prune_orphans()
    return graph


def build_readout_graph(conn, drugs: list[str], max_edges: int) -> Graph:
    """drug -> biomarker category, disease-relevant readouts only.

    PK is excluded: it is 76% of the annotated measures and says only that the
    drug was in the bloodstream, not that anything happened to the disease.
    """
    graph = Graph(max_edges)
    for drug, category, n_trials, n_measures in conn.execute(
        """
        SELECT drug, category, count(DISTINCT nct_id), count(*)
        FROM biomarker_measures
        WHERE drug IN ? AND is_disease_readout
        GROUP BY 1, 2
        ORDER BY 3 DESC
        """,
        [drugs],
    ).fetchall():
        drug_id = graph.node(f"drug:{drug}", "drug", drug, rank=drugs.index(drug) + 1)
        category_id = graph.node(f"readout:{category}", "readout", category)
        graph.edge(
            drug_id,
            category_id,
            "measured_by",
            f"{n_trials} trials",
            n_trials=n_trials,
            n_measures=n_measures,
        )
    graph.prune_orphans()
    return graph


def write_tables(graph: Graph, out_dir: Path, stem: str) -> None:
    node_fields = sorted({k for n in graph.nodes.values() for k in n})
    edge_fields = sorted({k for e in graph.edges for k in e})
    for name, fields, rows in (
        (f"{stem}_nodes.csv", node_fields, graph.nodes.values()),
        (f"{stem}_edges.csv", edge_fields, graph.edges),
    ):
        with (out_dir / name).open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


def write_graphml(graph: Graph, path: Path) -> None:
    node_keys = sorted({k for n in graph.nodes.values() for k in n} - {"id"})
    edge_keys = sorted({k for e in graph.edges for k in e} - {"source", "target"})

    def kind_of(value) -> str:
        return (
            "double" if isinstance(value, float) else "int" if isinstance(value, int) else "string"
        )

    node_types = dict.fromkeys(node_keys, "string")
    edge_types = dict.fromkeys(edge_keys, "string")
    for n in graph.nodes.values():
        for k, v in n.items():
            if k in node_types and v is not None:
                node_types[k] = kind_of(v)
    for e in graph.edges:
        for k, v in e.items():
            if k in edge_types and v is not None:
                edge_types[k] = kind_of(v)

    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
    ]
    for key in node_keys:
        out.append(
            f'  <key id="n_{key}" for="node" attr.name="{key}" attr.type="{node_types[key]}"/>'
        )
    for key in edge_keys:
        out.append(
            f'  <key id="e_{key}" for="edge" attr.name="{key}" attr.type="{edge_types[key]}"/>'
        )
    out.append('  <graph id="G" edgedefault="directed">')
    for node in graph.nodes.values():
        out.append(f'    <node id="{escape(node["id"])}">')
        for key in node_keys:
            value = node.get(key)
            if value is not None:
                out.append(f'      <data key="n_{key}">{escape(str(value))}</data>')
        out.append("    </node>")
    for i, edge in enumerate(graph.edges):
        out.append(
            f'    <edge id="e{i}" source="{escape(edge["source"])}" '
            f'target="{escape(edge["target"])}">'
        )
        for key in edge_keys:
            value = edge.get(key)
            if value is not None:
                out.append(f'      <data key="e_{key}">{escape(str(value))}</data>')
        out.append("    </edge>")
    out += ["  </graph>", "</graphml>"]
    path.write_text("\n".join(out) + "\n")


def _mermaid_label(text: str) -> str:
    return text.replace('"', "'").replace("\n", " ")


def write_mermaid(graph: Graph, path: Path, *, title: str) -> None:
    shape = {
        "drug": ('["', '"]'),
        "anchor_drug": ('["', '"]'),
        "target": ('(["', '"])'),
        "toxin": ('{{"', '"}}'),
        "pathway": ('["', '"]'),
        "disease": ('[["', '"]]'),
        "readout": ('("', '")'),
    }
    ids = {node_id: f"n{i}" for i, node_id in enumerate(graph.nodes)}
    out = [f"%% {title}", "flowchart LR"]
    for node_id, node in graph.nodes.items():
        open_, close = shape.get(node["kind"], ('["', '"]'))
        label = node["label"]
        if node["kind"] == "drug" and node.get("rank"):
            label = f"{node['rank']}. {label}"
        elif node["kind"] == "anchor_drug":
            label = f"{label} (approved for PD)"
        out.append(f"  {ids[node_id]}{open_}{_mermaid_label(label)}{close}")
    for edge in graph.edges:
        label = f"|{_mermaid_label(edge['label'])}|" if edge["label"] else ""
        out.append(f"  {ids[edge['source']]} -->{label} {ids[edge['target']]}")
    out += [
        "  classDef drug fill:#ddf4ff,stroke:#0969da,color:#0a3069;",
        "  classDef anchor fill:#eeeeee,stroke:#6e7781,color:#24292f,stroke-dasharray:4 3;",
        "  classDef target fill:#fff8c5,stroke:#9a6700,color:#4d2d00;",
        "  classDef toxin fill:#ffebe9,stroke:#cf222e,color:#5a0f14;",
        "  classDef pathway fill:#dafbe1,stroke:#1a7f37,color:#0a3622;",
        "  classDef disease fill:#fbefff,stroke:#8250df,color:#3b1a5c;",
        "  classDef readout fill:#f6f8fa,stroke:#6e7781,color:#24292f;",
    ]
    for kind, cls in (
        ("drug", "drug"),
        ("anchor_drug", "anchor"),
        ("target", "target"),
        ("toxin", "toxin"),
        ("pathway", "pathway"),
        ("disease", "disease"),
        ("readout", "readout"),
    ):
        members = [ids[i] for i, n in graph.nodes.items() if n["kind"] == kind]
        if members:
            out.append(f"  class {','.join(members)} {cls};")
    path.write_text("\n".join(out) + "\n")


def write_pathway_graph(conn, path: Path, max_edges: int) -> Graph:
    """Pathway-level convergence: toxin -> pathway <- enrichment, with drug counts."""
    graph = Graph(max_edges)
    disease = graph.node(f"disease:{PD_ROOT}", "disease", PD_ROOT)
    for (
        index,
        name,
        size,
        qvalue,
        fold,
        _n_toxins,
        support,
        toxins,
        targets,
        n_target_drugs,
    ) in conn.execute(
        """
        SELECT e.pathway_index, e.pathway_name, e.pathway_size, e.pathway_qvalue,
               e.fold_enrichment, e.n_toxins, e.toxin_support, e.toxins, e.toxin_targets,
               coalesce(tp.n_target_drugs, 0)
        FROM pathway_evidence e
        LEFT JOIN toxin_pathway tp USING (pathway_index)
        WHERE e.is_convergent
        ORDER BY e.toxin_support DESC
        """
    ).fetchall():
        pathway_id = graph.node(
            f"pathway:{name}",
            "pathway",
            name,
            pathway_index=index,
            pathway_size=size,
            qvalue=qvalue,
            fold_enrichment=fold,
            toxin_support=round(support or 0, 4),
            n_target_drugs=n_target_drugs,
        )
        graph.edge(pathway_id, disease, "enriched_in", f"q={qvalue:.2g}", qvalue=qvalue)
        for target in filter(None, (t.strip() for t in (targets or "").split(";"))):
            graph.edge(
                graph.node(f"target:{target}", "target", target),
                pathway_id,
                "member_of",
                "in pathway",
            )
        for toxin in filter(None, (t.strip() for t in (toxins or "").split(";"))):
            graph.edge(
                graph.node(f"toxin:{toxin}", "toxin", toxin), pathway_id, "implicates", "implicates"
            )
    graph.prune_orphans()
    write_mermaid(graph, path, title="Convergent PD pathways: toxin evidence + enrichment")
    return graph


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--database", default=str(DEFAULT_DB))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--max-edges", type=int, default=100)
    args = parser.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(args.database, read_only=True)

    ranked = conn.execute(RANKING_SQL, [args.top]).fetchall()
    columns = [d[0] for d in conn.description]
    ranking_path = out_dir / "pd_ranking_top10.csv"
    with ranking_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rank", *columns])
        for i, row in enumerate(ranked, start=1):
            writer.writerow([i, *row])
    drugs = [r[0] for r in ranked]

    anchors = [r[0] for r in conn.execute(ANCHOR_SQL).fetchall()]

    mechanism = build_mechanism_graph(conn, drugs, anchors, args.max_edges)
    write_tables(mechanism, out_dir, "pd_top10")
    write_graphml(mechanism, out_dir / "pd_top10.graphml")
    write_mermaid(
        mechanism,
        out_dir / "pd_top10_mechanism.mmd",
        title="Top PD repurposing candidates: toxin, target and pathway evidence",
    )

    readouts = build_readout_graph(conn, drugs, args.max_edges)
    write_mermaid(
        readouts,
        out_dir / "pd_top10_readouts.mmd",
        title="Disease-relevant biomarker readouts (PK excluded)",
    )

    pathways = write_pathway_graph(conn, out_dir / "pd_pathway_convergence.mmd", args.max_edges)
    conn.close()

    print(f"ranking: {len(ranked)} drugs -> {ranking_path}", file=sys.stderr)
    for name, graph in (("mechanism", mechanism), ("readouts", readouts), ("pathways", pathways)):
        note = f"  ({graph.dropped} edges dropped to the budget)" if graph.dropped else ""
        print(
            f"{name:10} {len(graph.nodes):3} nodes  {len(graph.edges):3} edges"
            f"  (cap {args.max_edges}){note}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
