"""Export the evidence behind one drug as a small, readable graph.

A ranking says a drug scored well. This says *why*, in the form a scientist can
disagree with:

    toxin  --implicates-->  target  <--acts on--  drug
                              |
                              +--member of-->  pathway  --enriched in-->  Parkinson's

Edges from the drug are coloured by direction, which is the claim most worth
arguing with: green where the drug opposes the toxin insult, red where it mimics
it, amber where the curation says the argument genuinely runs both ways. A red
edge is not a low score, it is a different hypothesis -- and seeing it is the
point of exporting the graph at all.

Three renderings, same content: JSON for machines, Graphviz DOT for publication
figures, Mermaid for anything that renders markdown.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

PD_ROOT = "Parkinson disease"

# direction -> (graphviz colour, mermaid class)
DIRECTION_STYLE = {
    "protective": ("#1a7f37", "protective"),
    "risk": ("#cf222e", "risk"),
    "ambiguous": ("#bf8700", "ambiguous"),
    "unknown": ("#6e7781", "unknown"),
}


@dataclass
class Node:
    id: str
    kind: str  # drug | target | toxin | pathway | disease | trial
    label: str
    attrs: dict[str, Any] = field(default_factory=dict)


@dataclass
class Edge:
    source: str
    target: str
    kind: str
    label: str = ""
    attrs: dict[str, Any] = field(default_factory=dict)


@dataclass
class EvidenceGraph:
    """The evidence for one drug, as nodes and edges."""

    drug: str
    nodes: list[Node]
    edges: list[Edge]
    summary: dict[str, Any]

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(
            {
                "drug": self.drug,
                "summary": self.summary,
                "nodes": [
                    {"id": n.id, "kind": n.kind, "label": n.label, **n.attrs} for n in self.nodes
                ],
                "edges": [
                    {
                        "source": e.source,
                        "target": e.target,
                        "kind": e.kind,
                        "label": e.label,
                        **e.attrs,
                    }
                    for e in self.edges
                ],
            },
            indent=indent,
        )

    def to_dot(self) -> str:
        shape = {
            "drug": ("box", "#ddf4ff"),
            "target": ("ellipse", "#fff8c5"),
            "toxin": ("octagon", "#ffebe9"),
            "pathway": ("box", "#dafbe1"),
            "disease": ("doubleoctagon", "#f6f8fa"),
            "trial": ("note", "#f6f8fa"),
        }
        out = [
            f'digraph "{_esc(self.drug)}" {{',
            "  rankdir=LR;",
            '  node [style="filled,rounded" fontname="Helvetica" fontsize=10];',
            '  edge [fontname="Helvetica" fontsize=9];',
        ]
        for n in self.nodes:
            form, fill = shape.get(n.kind, ("box", "#ffffff"))
            out.append(
                f'  "{_esc(n.id)}" [label="{_esc(n.label)}" shape={form} fillcolor="{fill}"];'
            )
        for e in self.edges:
            colour = DIRECTION_STYLE.get(e.attrs.get("direction", ""), ("#57606a", ""))[0]
            style = " style=dashed" if e.kind == "trial" else ""
            out.append(
                f'  "{_esc(e.source)}" -> "{_esc(e.target)}" '
                f'[label="{_esc(e.label)}" color="{colour}" fontcolor="{colour}"{style}];'
            )
        out.append("}")
        return "\n".join(out)

    def to_mermaid(self) -> str:
        ids = {n.id: f"n{i}" for i, n in enumerate(self.nodes)}
        shape = {
            "drug": ('["', '"]'),
            "target": ('(["', '"])'),
            "toxin": ('{{"', '"}}'),
            "pathway": ('["', '"]'),
            "disease": ('[["', '"]]'),
            "trial": ('("', '")'),
        }
        out = ["flowchart LR"]
        for n in self.nodes:
            open_, close = shape.get(n.kind, ('["', '"]'))
            out.append(f"  {ids[n.id]}{open_}{_mermaid_text(n.label)}{close}")
        for e in self.edges:
            arrow = "-.->" if e.kind == "trial" else "-->"
            label = f"|{_mermaid_text(e.label)}|" if e.label else ""
            out.append(f"  {ids[e.source]} {arrow}{label} {ids[e.target]}")
        # Colour the direction claim, which is the part worth arguing with.
        out += [
            "  classDef protective fill:#dafbe1,stroke:#1a7f37,color:#0a3622;",
            "  classDef risk fill:#ffebe9,stroke:#cf222e,color:#5a0f14;",
            "  classDef ambiguous fill:#fff8c5,stroke:#bf8700,color:#4d2d00;",
            "  classDef unknown fill:#f6f8fa,stroke:#6e7781,color:#24292f;",
        ]
        for direction, (_, cls) in DIRECTION_STYLE.items():
            members = [
                ids[e.target]
                for e in self.edges
                if e.kind == "acts_on" and e.attrs.get("direction") == direction
            ]
            if members:
                out.append(f"  class {','.join(sorted(set(members)))} {cls};")
        return "\n".join(out)


def _esc(text: str) -> str:
    return text.replace("\\", "\\\\").replace('"', '\\"')


def _mermaid_text(text: str) -> str:
    """Mermaid has no escape for a quote inside a quoted label; substitute it.

    Line breaks are ``\\n`` for Graphviz but ``<br/>`` for Mermaid, so the
    literal two-character sequence is translated rather than left to render raw.
    """
    return text.replace('"', "'").replace("\\n", "<br/>").replace("\n", " ")


def evidence_graph(
    conn,
    drug: str,
    *,
    enriched_only: bool = True,
    max_pathways: int = 6,
    include_trials: bool = True,
    max_trials: int = 5,
) -> EvidenceGraph:
    """Build the evidence graph for one drug.

    ``enriched_only`` keeps pathways at q < 0.05; without it the target's whole
    pathway membership comes along, which is usually noise. Drug name matching is
    case-insensitive but otherwise exact.
    """
    rows = conn.execute(
        """
        SELECT drug, gene, protein, action_type, direction, direction_confidence,
               direction_rationale, n_toxins, toxin_support, toxins,
               best_pathway_qvalue, n_enriched_pathways
        FROM drug_target_evidence
        WHERE lower(drug) = lower(?)
        ORDER BY gene
        """,
        [drug],
    ).fetchall()
    if not rows:
        raise LookupError(f"no drug named {drug!r} in drug_target_evidence")

    resolved = rows[0][0]
    nodes: list[Node] = [Node(f"drug:{resolved}", "drug", resolved)]
    edges: list[Edge] = []
    seen: set[str] = {f"drug:{resolved}"}

    def add(node: Node) -> str:
        if node.id not in seen:
            seen.add(node.id)
            nodes.append(node)
        return node.id

    genes = [r[1] for r in rows]
    for (
        _,
        gene,
        protein,
        action,
        direction,
        confidence,
        rationale,
        _n_toxins,
        support,
        toxins,
        best_q,
        n_enriched,
    ) in rows:
        gid = add(
            Node(
                f"target:{gene}",
                "target",
                gene,
                {
                    "protein": protein,
                    "direction": direction,
                    "direction_confidence": confidence,
                    "rationale": rationale,
                    "best_pathway_qvalue": best_q,
                    "n_enriched_pathways": n_enriched,
                    "toxin_support": support,
                },
            )
        )
        edges.append(
            Edge(
                f"drug:{resolved}",
                gid,
                "acts_on",
                (action or "mechanism unknown").lower(),
                {"direction": direction, "direction_confidence": confidence},
            )
        )
        for toxin in filter(None, (t.strip() for t in (toxins or "").split(";"))):
            tid = add(Node(f"toxin:{toxin}", "toxin", toxin))
            edges.append(Edge(tid, gid, "implicates", "toxin target"))

    # target -> pathway -> disease
    placeholders = ", ".join("?" for _ in genes)
    clause = "AND pathway_qvalue < 0.05" if enriched_only else ""
    pathways = conn.execute(
        f"""
        SELECT protein_name, pathway_name, pathway_index, min(pathway_qvalue) AS q
        FROM pd_tree
        WHERE protein_name IN ({placeholders}) {clause}
        GROUP BY 1, 2, 3
        ORDER BY q
        LIMIT {int(max_pathways)}
        """,
        genes,
    ).fetchall()
    for gene, pathway, _index, qvalue in pathways:
        pid = add(Node(f"pathway:{pathway}", "pathway", pathway, {"qvalue": qvalue}))
        edges.append(Edge(f"target:{gene}", pid, "member_of", "in pathway"))
        did = add(Node(f"disease:{PD_ROOT}", "disease", PD_ROOT))
        edges.append(Edge(pid, did, "enriched_in", f"q={qvalue:.2g}"))

    if include_trials:
        for nct, phase, status, conditions in conn.execute(
            """
            SELECT nct_id, phase, status, conditions
            FROM trials
            WHERE lower(drug) = lower(?) AND is_neuro
            ORDER BY is_parkinsons DESC, phase DESC
            LIMIT ?
            """,
            [resolved, int(max_trials)],
        ).fetchall():
            tid = add(
                Node(
                    f"trial:{nct}",
                    "trial",
                    f"{nct}\\n{conditions[:40]}",
                    {"phase": phase, "status": status, "conditions": conditions},
                )
            )
            edges.append(Edge(f"drug:{resolved}", tid, "trial", phase or "trial"))

    candidate = conn.execute(
        "SELECT verdict, priority_score, n_protective, n_risk, n_ambiguous, n_trials, "
        "n_pd_trials FROM drug_candidates WHERE lower(drug) = lower(?)",
        [resolved],
    ).fetchone()
    summary = {
        "targets": genes,
        "verdict": candidate[0] if candidate else None,
        "priority_score": candidate[1] if candidate else None,
        "n_protective": candidate[2] if candidate else None,
        "n_risk": candidate[3] if candidate else None,
        "n_ambiguous": candidate[4] if candidate else None,
        "n_trials": candidate[5] if candidate else None,
        "n_pd_trials": candidate[6] if candidate else None,
        "directions": {
            r[1]: {"action": r[3], "direction": r[4], "confidence": r[5], "rationale": r[6]}
            for r in rows
        },
    }
    return EvidenceGraph(resolved, nodes, edges, summary)
