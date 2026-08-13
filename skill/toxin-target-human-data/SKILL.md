---
name: toxin-target-human-data
description: >
  Mine the literature and chemical databases for environmental toxins associated with a
  disease, then enrich each toxin with human-relevant safety data and its biological
  (gene/protein) targets. Queries PubMed (NCBI E-utilities) for disease + environmental
  exposure papers and extracts a ranked toxin list; resolves each toxin in PubChem (PUG REST
  + PUG View) for human toxicity / clinical / epidemiology sections; pulls literature
  co-occurring gene targets (PubChem ChemicalGeneSymbolNeighbor) with gene-symbol acronyms;
  and merges everything into one integrated CSV. Worked example: Parkinson's disease
  environmental toxins. Disease-agnostic — swap the PubMed query and toxin keywords.
  Only needs Python 3 + requests.
---

# Toxin → Human Data → Target Mining

Turn a disease name into a ranked list of associated environmental toxins, each annotated
with **human-relevant evidence** and its **biological gene/protein targets**. Worked
example: Parkinson's disease (PD). The method is disease-agnostic — change the PubMed query
and the toxin keyword list.

## Pipeline

```
PubMed (E-utilities)          PubChem PUG REST/View            PubChem link_db
        │                              │                              │
  esearch + efetch            name → CID → full record        ChemicalGeneSymbolNeighbor
        │                              │                              │
  toxin list + counts   ──►    human-data sections     ──►     gene targets + acronyms
        └──────────────────────────────┴──────────────────────────────┘
                                        ▼
                        one integrated CSV (per toxin)
```

Three stages, each a script in `scripts/`:

| Stage | Script | Source | Output |
|-------|--------|--------|--------|
| 1. Toxin discovery | `pubmed_toxins.py` | PubMed E-utilities | papers CSV + toxin mention counts |
| 2. Human data | `pubchem_lookup.py` | PubChem PUG REST / PUG View | human toxicity / clinical / epi sections per compound |
| 3. Integrate + targets | `integrate_toxin_human_data.py` | ties stages 1–2 together + PubChem `link_db` | one merged CSV with gene targets |

## Core principle

A toxin appearing in the literature is a **hypothesis**, not a mechanism. Weight comes from
stacking independent signals:

- **Literature frequency** — how often the toxin co-occurs with the disease (PubMed).
- **Human relevance** — does PubChem carry human toxicity excerpts, exposure routes,
  clinical-trial or epidemiology sections (not just animal LD50s)?
- **Target plausibility** — do the toxin's literature gene targets include disease-relevant
  biology? (For PD, e.g. `SNCA`, `TH`, `SLC6A3` show up for MPTP and rotenone — a strong
  signal; a toxin whose only targets are generic stress enzymes is weaker.)

## How to run

Default worked example (Parkinson's disease):

```bash
cd scripts
# Stage 1: discover toxins from PubMed
python3 pubmed_toxins.py --out pd_toxin_papers.csv
# Stage 3 calls stage 2 internally; it reads stage 1's CSV and writes the merged file
python3 integrate_toxin_human_data.py --in pd_toxin_papers.csv --out pd_toxin_human_data_integrated.csv
```

A pre-built result for PD ships in `../data/pd_toxin_human_data_integrated.csv`.

Stage 2 can also be run standalone to inspect one compound:

```bash
python3 pubchem_lookup.py --name paraquat --full
```

## Adapting to another disease

1. In `pubmed_toxins.py`, edit `DEFAULT_QUERY` (swap the disease MeSH/title terms) and, if
   needed, extend `TOXIN_KEYWORDS` with candidate exposures for that disease.
2. Re-run the two commands above with new output names.
3. Nothing else changes — stages 2 and 3 are disease-agnostic (they operate on whatever
   compound names stage 1 surfaces).

## Output schema (integrated CSV)

| Column | Meaning |
|--------|---------|
| `Toxin` | toxin/compound name surfaced from PubMed abstracts |
| `PubMedMentions` | # of fetched PD papers mentioning it |
| `CID` | PubChem Compound ID (blank if the name is a class/exposure, not one molecule) |
| `HasHumanData` | True if PubChem has human toxicity/clinical/epi sections |
| `MatchedSections` | which human-data section headings matched |
| `NumGeneTargets` | # literature co-occurring gene/protein targets |
| `TopGeneTargets` | top targets as `full name [SYMBOL] (article_count)`; symbol is a gene acronym (e.g. `SNCA`) or an EC number for generic enzyme families |
| `SampleExcerpt` | one short human-data excerpt |

## Notes & limits

- **Rate limits.** PubMed allows 3 req/s without a key, 10 with `--api-key`. PubChem calls
  are spaced with short sleeps; keep them.
- **Unresolved names.** Exposure classes ("air pollution", "particulate matter",
  "organophosphate") don't map to a single PubChem CID and come back with blank target /
  human-data columns — expected, not an error.
- **Targets are literature co-occurrence**, not assay-confirmed binding. High article counts
  mean a well-studied association, which usually (not always) reflects human/experimental
  evidence. Treat as a lead-generation signal.
- CTD (Comparative Toxicogenomics Database) is the richer curated toxin→gene source but is
  now behind an anti-bot wall; this skill uses the PubChem-native `link_db` endpoint instead
  for reproducibility.
