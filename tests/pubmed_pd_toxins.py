#!/usr/bin/env python3
"""
pubmed_pd_toxins.py

Pull literature from the PubMed API (NCBI E-utilities,
https://www.ncbi.nlm.nih.gov/books/NBK25501/) on environmental toxins
associated with increased Parkinson's disease (PD) risk.

Uses:
  - esearch  -> get a list of PMIDs matching a query
  - efetch   -> get title/abstract/journal/date/authors for those PMIDs

Examples
--------
# Default query: environmental toxins + Parkinson's disease risk
python3 pubmed_pd_toxins.py --out pd_toxin_papers.csv

# Custom query, more results
python3 pubmed_pd_toxins.py --query "paraquat AND Parkinson disease" --max 200 --out paraquat.csv

# Add your own NCBI API key (raises rate limit from 3 to 10 req/sec)
python3 pubmed_pd_toxins.py --api-key YOUR_KEY --out pd_toxin_papers.csv
"""

import argparse
import csv
import sys
import time
import xml.etree.ElementTree as ET

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_QUERY = (
    '("Parkinson Disease"[MeSH Terms] OR "Parkinson\'s disease"[Title/Abstract]) '
    'AND ("environmental exposure"[MeSH Terms] OR pesticide[Title/Abstract] OR '
    "herbicide[Title/Abstract] OR solvent[Title/Abstract] OR toxin[Title/Abstract] "
    'OR toxicant[Title/Abstract] OR "heavy metal"[Title/Abstract]) '
    'AND ("risk"[Title/Abstract] OR "association"[Title/Abstract])'
)

# Known/candidate PD-linked environmental toxins to flag when scanning abstracts.
# Extend this list as needed.
TOXIN_KEYWORDS = [
    "paraquat",
    "rotenone",
    "chlorpyrifos",
    "organochlorine",
    "organophosphate",
    "trichloroethylene",
    "tce",
    "perchloroethylene",
    "tetrachloroethylene",
    "pce",
    "manganese",
    "lead",
    "mercury",
    "cadmium",
    "arsenic",
    "mptp",
    "agent orange",
    "dioxin",
    "pcb",
    "polychlorinated biphenyl",
    "air pollution",
    "particulate matter",
    "pm2.5",
    "traffic-related",
    "dieldrin",
    "ddt",
    "atrazine",
    "glyphosate",
    "permethrin",
]


def esearch(query, retmax=100, api_key=None):
    """Return a list of PMIDs matching the query."""
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": retmax,
        "retmode": "json",
        "sort": "relevance",
    }
    if api_key:
        params["api_key"] = api_key
    resp = requests.get(f"{EUTILS}/esearch.fcgi", params=params, timeout=20)
    resp.raise_for_status()
    return resp.json()["esearchresult"]["idlist"]


def efetch_batch(pmids, api_key=None):
    """Fetch title/abstract/journal/date/authors for a batch of PMIDs (XML)."""
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
    }
    if api_key:
        params["api_key"] = api_key
    resp = requests.get(f"{EUTILS}/efetch.fcgi", params=params, timeout=30)
    resp.raise_for_status()
    return ET.fromstring(resp.content)


def parse_article(article_elem):
    """Extract fields of interest from one <PubmedArticle> XML element."""
    pmid = article_elem.findtext(".//PMID", default="")
    title = article_elem.findtext(".//ArticleTitle", default="")
    journal = article_elem.findtext(".//Journal/Title", default="")
    year = article_elem.findtext(".//JournalIssue/PubDate/Year") or article_elem.findtext(
        ".//JournalIssue/PubDate/MedlineDate", default=""
    )
    abstract_parts = [(elem.text or "") for elem in article_elem.findall(".//AbstractText")]
    abstract = " ".join(abstract_parts).strip()

    authors = []
    for author in article_elem.findall(".//AuthorList/Author"):
        last = author.findtext("LastName")
        fore = author.findtext("ForeName")
        if last:
            authors.append(f"{fore} {last}".strip() if fore else last)
    authors_str = "; ".join(authors[:5]) + (" et al." if len(authors) > 5 else "")

    text_blob = f"{title} {abstract}".lower()
    toxins_found = sorted({kw for kw in TOXIN_KEYWORDS if kw in text_blob})

    return {
        "PMID": pmid,
        "Title": title,
        "Journal": journal,
        "Year": year,
        "Authors": authors_str,
        "ToxinsMentioned": "; ".join(toxins_found),
        "Abstract": abstract,
        "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    }


def fetch_pubmed(query, max_results=100, api_key=None, batch_size=50):
    pmids = esearch(query, retmax=max_results, api_key=api_key)
    print(f"Found {len(pmids)} matching PMIDs")

    rows = []
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        print(f"  fetching records {i + 1}-{i + len(batch)} ...")
        root = efetch_batch(batch, api_key=api_key)
        for article_elem in root.findall(".//PubmedArticle"):
            rows.append(parse_article(article_elem))
        time.sleep(0.34 if not api_key else 0.11)  # respect NCBI rate limits
    return rows


def summarize_toxins(rows):
    """Tally how often each toxin keyword appears across the fetched abstracts."""
    counts = {}
    for row in rows:
        for toxin in row["ToxinsMentioned"].split("; "):
            if toxin:
                counts[toxin] = counts.get(toxin, 0) + 1
    return sorted(counts.items(), key=lambda kv: kv[1], reverse=True)


def main():
    parser = argparse.ArgumentParser(
        description="Query PubMed for PD/environmental toxin literature"
    )
    parser.add_argument("--query", default=DEFAULT_QUERY, help="PubMed search query")
    parser.add_argument("--max", type=int, default=100, help="max number of records to fetch")
    parser.add_argument(
        "--api-key", default=None, help="NCBI API key (optional, raises rate limit)"
    )
    parser.add_argument("--out", default="pd_toxin_papers.csv", help="output CSV path")
    args = parser.parse_args()

    rows = fetch_pubmed(args.query, max_results=args.max, api_key=args.api_key)
    if not rows:
        print("No results found.")
        sys.exit(0)

    fieldnames = list(rows[0].keys())
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved {len(rows)} record(s) to {args.out}")

    print("\nToxin mention frequency across fetched abstracts:")
    for toxin, count in summarize_toxins(rows):
        print(f"  {toxin:25s} {count}")


if __name__ == "__main__":
    main()
