#!/usr/bin/env python3
"""
integrate_pd_toxin_human_data.py

Pipeline glue:
  1. Read the toxin list (+ PubMed mention counts) produced by
     pubmed_pd_toxins.py (default input: pd_toxin_papers.csv).
  2. For each toxin, run it through PubChem via pubchem_lookup.py to
     find human-relevant data (toxicity excerpts, clinical trials,
     epidemiology, human health effects, ...).
  3. Merge PubMed mention counts with PubChem human-data findings and
     export one integrated CSV.

Usage
-----
python3 integrate_pd_toxin_human_data.py \
    --in pd_toxin_papers.csv \
    --out pd_toxin_human_data_integrated.csv
"""

import argparse
import csv
import time
from collections import Counter

import requests

from pubchem_lookup import check_compound

LINK_DB = "https://pubchem.ncbi.nlm.nih.gov/link_db/link_db_server.cgi"


def gene_targets_for_cid(cid, top_n=8):
    """
    Fetch literature co-occurring gene/protein targets for a PubChem CID.

    Uses PubChem's ChemicalGeneSymbolNeighbor link set: genes/proteins that
    co-occur with the compound across the biomedical literature, ranked by
    article count (a proxy for how well the toxin->target link is supported
    by human/experimental studies). Returns (num_targets, top_targets_str).
    """
    if not cid:
        return 0, ""
    params = {
        "format": "JSON",
        "type": "ChemicalGeneSymbolNeighbor",
        "operation": "GetAllLinks",
        "id_1": str(cid),
    }
    try:
        resp = requests.get(LINK_DB, params=params, timeout=30,
                            headers={"User-Agent": "Mozilla/5.0"})
        resp.raise_for_status()
        links = resp.json().get("LinkDataSet", {}).get("LinkData", [])
    except (requests.RequestException, ValueError):
        return 0, ""

    targets = []
    for ld in links:
        ev = ld.get("Evidence", {}).get("ChemicalGeneSymbolNeighbor", {})
        name = ev.get("NeighborName", "").strip()
        if not name:
            continue
        # Gene symbol / acronym for easier search (e.g. SNCA, CAT, or an EC number)
        sym = ld.get("ID_2", {}).get("GeneSymbol", "").strip()
        if sym.lower().startswith("ec:"):
            sym = "EC:" + sym.split(":", 1)[1]
        elif sym:
            sym = sym.upper()
        targets.append((name, sym, ev.get("ArticleCount", 0)))
    targets.sort(key=lambda t: t[2], reverse=True)
    top = "; ".join(
        (f"{name} [{sym}] ({count})" if sym else f"{name} ({count})")
        for name, sym, count in targets[:top_n]
    )
    return len(targets), top


def toxin_counts_from_pubmed_csv(path):
    """Tally distinct toxins and their mention counts from the PubMed CSV."""
    counts = Counter()
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            for toxin in row.get("ToxinsMentioned", "").split("; "):
                toxin = toxin.strip()
                if toxin:
                    counts[toxin] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Integrate PubMed PD toxin list with PubChem human-data lookup"
    )
    parser.add_argument("--in", dest="infile", default="pd_toxin_papers.csv",
                        help="PubMed CSV from pubmed_pd_toxins.py")
    parser.add_argument("--out", default="pd_toxin_human_data_integrated.csv",
                        help="integrated output CSV path")
    args = parser.parse_args()

    counts = toxin_counts_from_pubmed_csv(args.infile)
    toxins = [t for t, _ in counts.most_common()]
    print(f"Loaded {len(toxins)} distinct toxins from {args.infile}\n")

    rows = []
    for toxin in toxins:
        print(f"Checking PubChem for '{toxin}' ...")
        result = check_compound(toxin, id_type="name")
        if result is None:
            # Name could not be resolved to a PubChem CID (e.g. "air pollution")
            result = {
                "Query": toxin, "CID": "", "HasHumanData": False,
                "MatchedSections": "", "SampleExcerpt": "",
            }
        n_targets, top_targets = gene_targets_for_cid(result["CID"])
        merged = {
            "Toxin": toxin,
            "PubMedMentions": counts[toxin],
            "CID": result["CID"],
            "HasHumanData": result["HasHumanData"],
            "MatchedSections": result["MatchedSections"],
            "NumGeneTargets": n_targets,
            "TopGeneTargets": top_targets,
            "SampleExcerpt": result["SampleExcerpt"],
        }
        rows.append(merged)
        status = "HAS human data" if result["HasHumanData"] else "no human data / unresolved"
        print(f"  CID {result['CID'] or '-'}: {status}; {n_targets} gene target(s)")
        time.sleep(0.3)  # be polite to the PubChem API

    fieldnames = ["Toxin", "PubMedMentions", "CID", "HasHumanData",
                  "MatchedSections", "NumGeneTargets", "TopGeneTargets",
                  "SampleExcerpt"]
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    n_human = sum(1 for r in rows if r["HasHumanData"])
    print(f"\nSaved {len(rows)} toxins to {args.out}")
    print(f"{n_human}/{len(rows)} toxins have human-relevant PubChem data.")


if __name__ == "__main__":
    main()
