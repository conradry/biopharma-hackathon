"""Map PD-associated environmental toxins onto PrimeKG proteins and pathways.

Step 1 -- parse ``pd_toxin_human_data_integrated.csv``, resolve its gene targets
to PrimeKG ``gene/protein`` nodes, and emit one row per (toxin, gene).

Step 2 -- weight each (toxin, gene) pair by specificity across toxins, then roll
the weights up to the pathways of ``pd_tree.parquet``.

The source ``TopGeneTargets`` field is a ``;``-separated list of
``name [TAG] (count)``, where TAG is either a gene symbol or an EC number, and
the count is a PubChem literature co-mention -- not an affinity. Three things
about it drive the design:

* The EC entries carry the mechanism. MPTP's causal target is monoamine oxidase,
  which appears only as ``EC:1.4.3.4``; symbol-only parsing loses it. EC_MAP
  below is hand-curated, and every expansion records how confident it is.
* The gene lists are cross-species (influenza PB2, bacterial MERA, rice
  LOC4326471). Nothing is trusted that fails to resolve against PrimeKG.
* A toxin often lists the same protein twice, once by symbol and once by EC
  (lead has both ``catalase [CAT]`` and ``catalase [EC:1.11.1.6]``). Rows are
  therefore aggregated to the (toxin, gene) grain; summing the raw list would
  double-count.

Usage:
    python scripts/build_toxin_targets.py [TOXIN_CSV] [KG_CSV] [OUT_DIR]
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from math import log
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from extract_pd_subgraph import DEFAULT_KG

DEFAULT_TOXINS = Path("data/pd_toxin_human_data_integrated.csv")
DEFAULT_TREE = Path("data/pd_tree.parquet")
DEFAULT_OUT = Path("data")

ENTRY = re.compile(r"^(?P<name>.+?)\s*\[(?P<tag>[^\]]+)\]\s*\((?P<count>\d+)\)$")

# EC class -> (human gene symbols, confidence). Confidence is recorded on every
# row rather than baked into a keep/drop decision, so a query can tighten the
# mapping without a rebuild:
#
#   exact       the EC names exactly one human gene
#   family      a small, fully enumerated human family (all members emitted)
#   complex     subunits of one complex; the toxin acts on the assembly
#   ambiguous   a large cross-species family -- no defensible human mapping
#   non_human   plant, fungal, bacterial or viral only
EC_MAP: dict[str, tuple[tuple[str, ...], str]] = {
    "EC:1.11.1.6": (("CAT",), "exact"),
    "EC:1.14.16.2": (("TH",), "exact"),
    "EC:3.1.1.7": (("ACHE",), "exact"),
    "EC:3.1.1.8": (("BCHE",), "exact"),
    "EC:4.2.1.24": (("ALAD",), "exact"),
    "EC:1.4.3.4": (("MAOA", "MAOB"), "family"),
    "EC:1.15.1.1": (("SOD1", "SOD2", "SOD3"), "family"),
    "EC:1.11.1.9": (("GPX1", "GPX2", "GPX3", "GPX4"), "family"),
    "EC:2.3.2.2": (("GGT1", "GGT5"), "family"),
    # Rotenone binds the ND1/NDUFS2 interface of complex I; the core catalytic
    # subunits stand in for the assembly. EC 1.6.99.3 itself is deprecated.
    "EC:1.6.99.3": (
        ("NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2"),
        "complex",
    ),
    "EC:1.11.1.7": ((), "ambiguous"),  # peroxidase: MPO/LPO/EPX/TPO/PXDN + non-human
    "EC:1.14.14.1": ((), "ambiguous"),  # unspecific monooxygenase: the whole CYP set
    "EC:1.10.3.2": ((), "non_human"),  # laccase (fungal)
    "EC:1.10.3.9": ((), "non_human"),  # photosystem II (plant)
    "EC:1.11.1.13": ((), "non_human"),  # manganese peroxidase (fungal)
    "EC:1.11.1.14": ((), "non_human"),  # lignin peroxidase (fungal)
    "EC:1.12.7.2": ((), "non_human"),  # ferredoxin hydrogenase (bacterial)
    "EC:1.21.99.5": ((), "non_human"),  # tetrachloroethene dehalogenase (bacterial)
}

# Symbols the source spells in a way PrimeKG will not match. Kept explicit and
# tiny -- this is a correction list, not a synonym resolver.
SYMBOL_FIXES = {"FER1HCH": "FTH1"}

TARGET_SCHEMA = pa.schema(
    [
        ("toxin", pa.string()),
        ("cid", pa.int32()),
        ("pubmed_mentions", pa.int32()),
        ("gene_symbol", pa.string()),
        ("protein_index", pa.int32()),
        ("evidence_count", pa.float64()),
        ("n_source_entries", pa.int32()),
        ("source_entries", pa.string()),
        ("mapping_confidence", pa.string()),
        ("toxin_mapped_fraction", pa.float64()),
        ("in_pd_tree", pa.bool_()),
        ("is_pd_seed", pa.bool_()),
        ("tf", pa.float64()),
        ("idf", pa.float64()),
        ("specificity_weight", pa.float64()),
    ]
)

PATHWAY_SCHEMA = pa.schema(
    [
        ("pathway_index", pa.int32()),
        ("pathway_name", pa.string()),
        ("pathway_size", pa.int32()),
        ("pathway_qvalue", pa.float64()),
        ("pathway_fold_enrichment", pa.float64()),
        ("n_toxins", pa.int32()),
        ("n_toxin_targets", pa.int32()),
        ("toxin_support", pa.float64()),
        ("toxins", pa.string()),
        ("toxin_targets", pa.string()),
        ("n_drugs", pa.int32()),
        ("n_target_drugs", pa.int32()),
    ]
)


@dataclass
class Entry:
    """One parsed ``name [TAG] (count)`` item, after EC expansion."""

    toxin: str
    label: str
    tag: str
    symbol: str
    count: int
    confidence: str
    # How many genes this entry's tag expanded to. Its evidence is split that
    # many ways, so a family tag cannot outvote a specific symbol.
    expansion: int


@dataclass
class Toxin:
    name: str
    cid: int | None
    pubmed_mentions: int
    entries: list[Entry] = field(default_factory=list)
    unmapped: list[tuple[str, str]] = field(default_factory=list)
    # Raw co-mention total across every listed entry, mapped or not. The
    # denominator for how much of a toxin's evidence survived resolution.
    listed_count: int = 0


def parse_toxins(path: Path) -> list[Toxin]:
    """Parse the CSV and expand EC tags into human gene symbols."""
    toxins: list[Toxin] = []
    for row in csv.DictReader(path.open(newline="")):
        if not row["TopGeneTargets"]:
            continue
        toxin = Toxin(
            name=row["Toxin"],
            cid=int(row["CID"]) if row["CID"] else None,
            pubmed_mentions=int(row["PubMedMentions"]),
        )
        for chunk in row["TopGeneTargets"].split(";"):
            chunk = chunk.strip()
            match = ENTRY.match(chunk)
            if not match:
                toxin.unmapped.append((chunk, "unparsed"))
                continue
            tag, count = match["tag"], int(match["count"])
            toxin.listed_count += count
            if tag.startswith("EC:"):
                symbols, confidence = EC_MAP.get(tag, ((), "ambiguous"))
                if not symbols:
                    toxin.unmapped.append((chunk, confidence))
                    continue
            else:
                symbols, confidence = (SYMBOL_FIXES.get(tag, tag),), "exact"
            for symbol in symbols:
                toxin.entries.append(
                    Entry(
                        toxin.name,
                        match["name"],
                        tag,
                        symbol,
                        count,
                        confidence,
                        len(symbols),
                    )
                )
        toxins.append(toxin)
    return toxins


def resolve_symbols(kg_path: Path, symbols: set[str]) -> dict[str, int]:
    """One streaming pass for the PrimeKG index of each gene/protein symbol.

    Off-tree symbols (complex I subunits, say) are not in pd_tree.parquet, so the
    lookup has to go back to the source edge list.
    """
    found: dict[str, set[int]] = defaultdict(set)
    with kg_path.open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            if row[4] == "gene/protein" and row[5] in symbols:
                found[row[5]].add(int(row[2]))
            if row[9] == "gene/protein" and row[10] in symbols:
                found[row[10]].add(int(row[7]))
    for symbol, indices in sorted(found.items()):
        if len(indices) > 1:
            print(f"warning: {symbol} maps to {sorted(indices)}", file=sys.stderr)
    return {symbol: min(indices) for symbol, indices in found.items()}


def main() -> None:
    argv = sys.argv[1:]
    toxin_csv = Path(argv[0]) if len(argv) > 0 else DEFAULT_TOXINS
    kg_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_KG
    out_dir = Path(argv[2]) if len(argv) > 2 else DEFAULT_OUT

    toxins = parse_toxins(toxin_csv)
    dropped: dict[str, int] = defaultdict(int)
    for toxin in toxins:
        for _, reason in toxin.unmapped:
            dropped[reason] += 1
    print(
        f"toxins with targets: {len(toxins)}  "
        f"expanded entries: {sum(len(t.entries) for t in toxins)}  "
        f"dropped: {dict(dropped)}",
        file=sys.stderr,
    )

    symbols = {e.symbol for t in toxins for e in t.entries}
    index_of = resolve_symbols(kg_path, symbols)
    unresolved = sorted(symbols - index_of.keys())
    print(
        f"symbols: {len(symbols)}  resolved in PrimeKG: {len(index_of)}  "
        f"unresolved (assumed non-human): {unresolved}",
        file=sys.stderr,
    )

    tree = pq.read_table(DEFAULT_TREE).to_pylist()
    in_tree = {r["protein_index"] for r in tree}
    seeds = {r["protein_index"] for r in tree if r["protein_is_pd_seed"]}

    # Collapse to the (toxin, gene) grain. A toxin listing a protein by both
    # symbol and EC contributes both counts to one row.
    merged: dict[tuple[str, str], dict[str, object]] = {}
    for toxin in toxins:
        for entry in toxin.entries:
            if entry.symbol not in index_of:
                continue
            key = (toxin.name, entry.symbol)
            cell = merged.setdefault(key, {"count": 0.0, "labels": [], "confidences": set()})
            cell["count"] += entry.count / entry.expansion
            cell["labels"].append(f"{entry.label} [{entry.tag}] ({entry.count})")
            cell["confidences"].add(entry.confidence)

    # Specificity: share of a toxin's mapped evidence, damped by how many toxins
    # name the same gene. CAT and SOD appear for 6 of 10 toxins and are the
    # background of the tox literature, not a toxin-specific signal.
    n_toxins = len({t for t, _ in merged})
    doc_freq: dict[str, int] = defaultdict(int)
    for _, symbol in merged:
        doc_freq[symbol] += 1
    totals: dict[str, float] = defaultdict(float)
    for (toxin_name, _), cell in merged.items():
        totals[toxin_name] += cell["count"]

    by_toxin = {t.name: t for t in toxins}
    rows = []
    for (toxin_name, symbol), cell in merged.items():
        protein_index = index_of[symbol]
        tf = cell["count"] / totals[toxin_name]
        idf = log(1 + n_toxins / doc_freq[symbol])
        order = ("exact", "complex", "family")
        confidence = min(cell["confidences"], key=lambda c: order.index(c))
        rows.append(
            {
                "toxin": toxin_name,
                "cid": by_toxin[toxin_name].cid,
                "pubmed_mentions": by_toxin[toxin_name].pubmed_mentions,
                "gene_symbol": symbol,
                "protein_index": protein_index,
                "evidence_count": round(cell["count"], 4),
                "n_source_entries": len(cell["labels"]),
                "source_entries": "; ".join(sorted(cell["labels"])),
                "mapping_confidence": confidence,
                # tf is a share of *mapped* evidence. Where most of a toxin's
                # list was non-human, the survivors absorb the whole share, so
                # a low fraction here means an inflated tf.
                "toxin_mapped_fraction": round(
                    totals[toxin_name] / by_toxin[toxin_name].listed_count, 4
                ),
                "in_pd_tree": protein_index in in_tree,
                "is_pd_seed": protein_index in seeds,
                "tf": round(tf, 6),
                "idf": round(idf, 6),
                "specificity_weight": round(tf * idf, 6),
            }
        )
    rows.sort(key=lambda r: (r["toxin"], -r["specificity_weight"]))

    out_dir.mkdir(parents=True, exist_ok=True)
    target_path = out_dir / "pd_toxin_target.parquet"
    pq.write_table(
        pa.table(
            {f.name: pa.array([r[f.name] for r in rows], f.type) for f in TARGET_SCHEMA},
            schema=TARGET_SCHEMA,
        ),
        target_path,
        compression="zstd",
    )

    # Roll the weights up to pathways. Only proteins already in the tree can
    # contribute; off-tree targets are reported but have no pathway to land on.
    weight_of: dict[int, list[tuple[str, float]]] = defaultdict(list)
    for row in rows:
        if row["in_pd_tree"]:
            weight_of[row["protein_index"]].append((row["toxin"], row["specificity_weight"]))

    pathways: dict[int, dict[str, object]] = {}
    for r in tree:
        p = pathways.setdefault(
            r["pathway_index"],
            {
                "pathway_index": r["pathway_index"],
                "pathway_name": r["pathway_name"],
                "pathway_size": r["pathway_size"],
                "pathway_qvalue": r["pathway_qvalue"],
                "pathway_fold_enrichment": r["pathway_fold_enrichment"],
                "proteins": set(),
                "drugs": set(),
                "target_drugs": set(),
            },
        )
        p["proteins"].add((r["protein_index"], r["protein_name"]))
        if r["drug_index"] is not None:
            p["drugs"].add(r["drug_index"])
            if r["drug_relation"] == "target":
                p["target_drugs"].add(r["drug_index"])

    pathway_rows = []
    for p in pathways.values():
        hits = [(idx, name) for idx, name in p["proteins"] if idx in weight_of]
        if not hits:
            continue
        toxin_names = {t for idx, _ in hits for t, _ in weight_of[idx]}
        support = sum(w for idx, _ in hits for _, w in weight_of[idx])
        pathway_rows.append(
            {
                "pathway_index": p["pathway_index"],
                "pathway_name": p["pathway_name"],
                "pathway_size": p["pathway_size"],
                "pathway_qvalue": p["pathway_qvalue"],
                "pathway_fold_enrichment": p["pathway_fold_enrichment"],
                "n_toxins": len(toxin_names),
                "n_toxin_targets": len(hits),
                "toxin_support": round(support, 6),
                "toxins": "; ".join(sorted(toxin_names)),
                "toxin_targets": "; ".join(sorted(name for _, name in hits)),
                "n_drugs": len(p["drugs"]),
                "n_target_drugs": len(p["target_drugs"]),
            }
        )
    pathway_rows.sort(key=lambda r: (-r["toxin_support"], r["pathway_qvalue"]))

    pathway_path = out_dir / "pd_toxin_pathway.parquet"
    pq.write_table(
        pa.table(
            {f.name: pa.array([r[f.name] for r in pathway_rows], f.type) for f in PATHWAY_SCHEMA},
            schema=PATHWAY_SCHEMA,
        ),
        pathway_path,
        compression="zstd",
    )

    off_tree = sorted({r["gene_symbol"] for r in rows if not r["in_pd_tree"]})
    significant = sum(1 for r in pathway_rows if r["pathway_qvalue"] < 0.05)
    print(
        f"(toxin, gene) rows: {len(rows)}  genes: {len(doc_freq)}  "
        f"off-tree genes: {len(off_tree)} {off_tree} -> {target_path}",
        file=sys.stderr,
    )
    print(
        f"pathways with toxin support: {len(pathway_rows)} "
        f"({significant} at q<0.05) -> {pathway_path}",
        file=sys.stderr,
    )
    # toxin_support is a raw sum, so it rewards large pathways the same way a
    # drug-count ranking rewards promiscuity. Enrichment controls for size, so
    # the convergent shortlist is q<0.05 first, support second.
    print("\nconvergent pathways (q<0.05, ranked by toxin support):", file=sys.stderr)
    for r in (p for p in pathway_rows if p["pathway_qvalue"] < 0.05):
        print(
            f"  {r['toxin_support']:7.3f}  q={r['pathway_qvalue']:<8.3g} "
            f"{r['n_toxins']} toxins  {r['n_target_drugs']:4} target drugs  "
            f"{r['pathway_name']}  [{r['toxin_targets']}]",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
