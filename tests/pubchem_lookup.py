#!/usr/bin/env python3
"""
pubchem_human_data.py

For a given compound, pull PubChem data and check whether PubChem's
compound record contains human-relevant data (human toxicity excerpts,
clinical trials, epidemiology, human health effects, etc.).

Uses two PubChem endpoints:
  - PUG REST  (https://pubchem.ncbi.nlm.nih.gov/rest/pug)      -> name -> CID, basic properties
  - PUG View  (https://pubchem.ncbi.nlm.nih.gov/rest/pug_view) -> full annotated record
                                                                   (toxicity, safety, literature, etc.)

Examples
--------
python3 pubchem_human_data.py --name paraquat
python3 pubchem_human_data.py --cid 15939
python3 pubchem_human_data.py --names paraquat rotenone manganese --out human_data_summary.csv
"""

import argparse
import csv
import sys
import time

import requests

PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

# Section headings (case-insensitive substring match) that indicate
# human-relevant data when present in a PubChem compound record.
HUMAN_DATA_HEADINGS = [
    "human toxicity excerpts",
    "human health effects",
    "clinical trials",
    "epidemiology",
    "exposure routes",
    "populations at special risk",
    "drug and medication information",
    "human metabolite information",
    "adverse effects",
]


def name_to_cid(name):
    """Resolve a compound name to a PubChem CID using PUG REST."""
    url = f"{PUG_REST}/compound/name/{name}/cids/JSON"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.json()["IdentifierList"]["CID"][0]


def get_full_record(cid, retries=3, pause=0.5):
    """Fetch the full annotated PUG View record for a CID."""
    url = f"{PUG_VIEW}/data/compound/{cid}/JSON"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as err:
            last_err = err
            time.sleep(pause * attempt)
    print(f"  ! failed to fetch full record for CID {cid}: {last_err}", file=sys.stderr)
    return None


def extract_strings(information_list, limit=2):
    """Pull a few short text excerpts out of a Section's Information list."""
    excerpts = []
    for info in information_list or []:
        value = info.get("Value", {})
        for item in value.get("StringWithMarkup", []):
            text = item.get("String", "").strip()
            if text:
                excerpts.append(text)
            if len(excerpts) >= limit:
                return excerpts
    return excerpts


def walk_sections(sections, path=None):
    """
    Recursively walk a PUG View Section tree, yielding
    (heading_path, information_list) for every section that has an
    Information list attached.
    """
    path = path or []
    for sec in sections or []:
        heading = sec.get("TOCHeading", "")
        new_path = [*path, heading]
        if sec.get("Information"):
            yield new_path, sec["Information"]
        if sec.get("Section"):
            yield from walk_sections(sec["Section"], new_path)


def find_human_data(record):
    """
    Scan a full PUG View record for sections matching HUMAN_DATA_HEADINGS.
    Returns a dict: heading -> list of short excerpts.
    """
    if not record:
        return {}
    root_sections = record.get("Record", {}).get("Section", [])
    matches = {}
    for path, info_list in walk_sections(root_sections):
        heading = path[-1].lower()
        for keyword in HUMAN_DATA_HEADINGS:
            if keyword in heading:
                excerpts = extract_strings(info_list)
                if excerpts:
                    matches[" > ".join(path)] = excerpts
                break
    return matches


def check_compound(identifier, id_type="name"):
    """Full pipeline for one compound: resolve CID, fetch record, scan for human data."""
    if id_type == "name":
        cid = name_to_cid(identifier)
        if cid is None:
            print(f"  ! could not resolve CID for '{identifier}'")
            return None
    else:
        cid = identifier

    record = get_full_record(cid)
    matches = find_human_data(record)

    return {
        "Query": identifier,
        "CID": cid,
        "HasHumanData": bool(matches),
        "MatchedSections": "; ".join(matches.keys()) if matches else "",
        "SampleExcerpt": (next(iter(matches.values()))[0][:300] if matches else ""),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Check PubChem for human-relevant data (toxicity, clinical trials, epidemiology) on a compound"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--name", help="single compound name")
    group.add_argument("--cid", help="single PubChem CID")
    group.add_argument("--names", nargs="+", help="multiple compound names")
    parser.add_argument("--out", help="write summary to this CSV file")
    parser.add_argument(
        "--full",
        action="store_true",
        help="print all matched excerpts in full, not just a sample",
    )
    args = parser.parse_args()

    if args.name:
        queries, id_type = [args.name], "name"
    elif args.cid:
        queries, id_type = [args.cid], "cid"
    else:
        queries, id_type = args.names, "name"

    results = []
    for q in queries:
        print(f"Checking {q} ...")
        result = check_compound(q, id_type=id_type)
        if result:
            results.append(result)
            status = "HAS human data" if result["HasHumanData"] else "no human data sections found"
            print(f"  CID {result['CID']}: {status}")
            if result["MatchedSections"]:
                print(f"  Sections: {result['MatchedSections']}")
        time.sleep(0.3)  # be polite to the API

    if not results:
        print("No results.")
        return

    if args.out:
        fieldnames = list(results[0].keys())
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSaved summary to {args.out}")

    if args.full:
        for q in queries:
            cid = name_to_cid(q) if id_type == "name" else q
            record = get_full_record(cid)
            matches = find_human_data(record)
            if matches:
                print(f"\n=== {q} (CID {cid}) — full excerpts ===")
                for heading, excerpts in matches.items():
                    print(f"\n[{heading}]")
                    for ex in excerpts:
                        print(f"  - {ex}")


if __name__ == "__main__":
    main()
