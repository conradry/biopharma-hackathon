"""Command line entry points for the PD knowledge base.

pdkb-build                          build pd_kb.duckdb
pdkb-rank --limit 20                the prioritised candidate list
pdkb-graph Rasagiline --format dot  why one drug looks promising
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from biopharma_hackathon.pdkb.build import DEFAULT_DB, connect
from biopharma_hackathon.pdkb.graph import evidence_graph


def _table(rows: list[tuple], headers: tuple[str, ...]) -> str:
    cells = [[f"{v:.3f}" if isinstance(v, float) else str(v) for v in r] for r in rows]
    widths = [
        max(len(h), *(len(c[i]) for c in cells)) if cells else len(h) for i, h in enumerate(headers)
    ]
    line = "  ".join(h.ljust(w) for h, w in zip(headers, widths, strict=True))
    out = [line, "  ".join("-" * w for w in widths)]
    out += ["  ".join(c.ljust(w) for c, w in zip(row, widths, strict=True)) for row in cells]
    return "\n".join(out)


def rank(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Rank repurposing candidates.")
    parser.add_argument("--database", default=str(DEFAULT_DB))
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument(
        "--verdict",
        help="Filter to one verdict, e.g. candidate, mechanistic_risk, "
        "established_pd_therapy, direction_ambiguous, direction_unknown.",
    )
    parser.add_argument(
        "--novel-only",
        action="store_true",
        help="Drop drugs that already have Parkinson's trials.",
    )
    args = parser.parse_args(argv)

    where = []
    params: list[object] = []
    if args.verdict:
        where.append("verdict = ?")
        params.append(args.verdict)
    if args.novel_only:
        where.append("n_pd_trials = 0")
    clause = f"WHERE {' AND '.join(where)}" if where else ""

    with connect(args.database) as conn:
        rows = conn.execute(
            f"""
            SELECT drug, targets, verdict, priority_score, direction_component,
                   pathway_component, toxin_component, n_trials, n_pd_trials
            FROM drug_candidates {clause}
            ORDER BY priority_score DESC LIMIT ?
            """,
            [*params, args.limit],
        ).fetchall()
    print(
        _table(
            rows,
            (
                "drug",
                "targets",
                "verdict",
                "score",
                "direction",
                "pathway",
                "toxin",
                "trials",
                "pd_trials",
            ),
        )
    )
    return 0


def graph(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export one drug's evidence graph.")
    parser.add_argument("drug")
    parser.add_argument("--database", default=str(DEFAULT_DB))
    parser.add_argument("--format", choices=("json", "dot", "mermaid"), default="mermaid")
    parser.add_argument("--out", help="Write to this path instead of stdout.")
    parser.add_argument(
        "--all-pathways",
        action="store_true",
        help="Include pathways that are not significantly enriched.",
    )
    parser.add_argument("--max-pathways", type=int, default=6)
    parser.add_argument("--no-trials", action="store_true")
    args = parser.parse_args(argv)

    with connect(args.database) as conn:
        try:
            evidence = evidence_graph(
                conn,
                args.drug,
                enriched_only=not args.all_pathways,
                max_pathways=args.max_pathways,
                include_trials=not args.no_trials,
            )
        except LookupError as exc:
            print(exc, file=sys.stderr)
            return 1

    text = {"json": evidence.to_json, "dot": evidence.to_dot, "mermaid": evidence.to_mermaid}[
        args.format
    ]()
    if args.out:
        Path(args.out).write_text(text + "\n")
        print(f"-> {args.out}", file=sys.stderr)
    else:
        print(text)
    return 0
