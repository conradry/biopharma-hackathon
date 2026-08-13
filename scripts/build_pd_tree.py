"""Build the Parkinson's disease repurposing tree from PrimeKG.

    Parkinson disease            (root)
      └── pathway                (283)
            └── protein          (target; the pathway's members)
                  └── drug       (drug_protein edge)

PrimeKG has no drug_pathway relation, so the protein layer is the only join
between drugs and pathways -- a drug hangs off the tree once per (pathway,
protein) pair it reaches.

Emitted as one row per root-to-leaf path, denormalized: every row carries the
whole ancestry, so any level rolls up with a GROUP BY. Nothing is pre-filtered;
the columns that would narrow the tree (``pathway_qvalue``,
``protein_is_pd_seed``, ``drug_relation``) ride along so you can cut it in the
query instead.

Proteins with no drug still get a row, with the drug columns null.

Usage:
    python scripts/build_pd_tree.py [KG_CSV] [OUT_PARQUET]
"""

from __future__ import annotations

import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from extract_pd_subgraph import DEFAULT_KG, PD_INDEX, resolve_nodes

DEFAULT_OUT = Path("data/pd_tree.parquet")
PD_NAME = "Parkinson disease"

SCHEMA = pa.schema(
    [
        ("disease_index", pa.int32()),
        ("disease_name", pa.string()),
        ("pathway_index", pa.int32()),
        ("pathway_name", pa.string()),
        ("pathway_size", pa.int32()),
        ("pathway_seed_overlap", pa.int32()),
        ("pathway_fold_enrichment", pa.float64()),
        ("pathway_qvalue", pa.float64()),
        ("protein_index", pa.int32()),
        ("protein_name", pa.string()),
        ("protein_is_pd_seed", pa.bool_()),
        ("drug_index", pa.int32()),
        ("drug_id", pa.string()),
        ("drug_name", pa.string()),
        ("drug_relation", pa.string()),
        ("depth", pa.int32()),
    ]
)


def main() -> None:
    kg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_KG
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUT

    scan = resolve_nodes(kg_path, with_drugs=True)
    rows: list[dict[str, object]] = []

    for pathway in sorted(scan.pathways, key=int):
        stats = scan.stats[pathway]
        branch = {
            "disease_index": int(PD_INDEX),
            "disease_name": PD_NAME,
            "pathway_index": int(pathway),
            "pathway_name": stats["pathway_name"],
            "pathway_size": stats["pathway_size"],
            "pathway_seed_overlap": stats["seed_overlap"],
            "pathway_fold_enrichment": stats["fold_enrichment"] or None,
            "pathway_qvalue": stats["qvalue"],
        }
        for protein in sorted(scan.pathway_members[pathway], key=int):
            leaf = {
                **branch,
                "protein_index": int(protein),
                "protein_name": scan.protein_names.get(protein),
                "protein_is_pd_seed": protein in scan.seed_proteins,
            }
            drugs = scan.drugs_by_protein.get(protein, {})
            if not drugs:
                rows.append(
                    {
                        **leaf,
                        "drug_index": None,
                        "drug_id": None,
                        "drug_name": None,
                        "drug_relation": None,
                        "depth": 3,
                    }
                )
                continue
            for (drug, display), (drug_id, drug_name) in sorted(
                drugs.items(), key=lambda kv: (int(kv[0][0]), kv[0][1])
            ):
                rows.append(
                    {
                        **leaf,
                        "drug_index": int(drug),
                        "drug_id": drug_id,
                        "drug_name": drug_name,
                        "drug_relation": display,
                        "depth": 4,
                    }
                )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.table(
        {f.name: pa.array([r[f.name] for r in rows], type=f.type) for f in SCHEMA},
        schema=SCHEMA,
    )
    pq.write_table(table, out_path, compression="zstd")

    drugs = {r["drug_index"] for r in rows if r["drug_index"] is not None}
    targets = {r["drug_index"] for r in rows if r["drug_relation"] == "target"}
    print(
        f"rows: {len(rows)}  pathways: {len(scan.pathways)}  "
        f"proteins: {len(scan.proteins)}  drugs: {len(drugs)} "
        f"({len(targets)} via a 'target' edge)  -> {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
