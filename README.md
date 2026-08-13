# biopharma-hackathon

## Setup

Requires [uv](https://docs.astral.sh/uv/).

```bash
uv sync                    # create .venv and install deps (incl. dev group)
uv run pre-commit install  # enable git hooks
```

## Common commands

```bash
uv run pytest              # run tests
uv run ruff check --fix .  # lint
uv run ruff format .       # format
uv add <pkg>               # add a runtime dependency
uv add --dev <pkg>         # add a dev dependency
uv run pre-commit run --all-files
```

Package code lives in `src/biopharma_hackathon/`, tests in `tests/`.

## Published data

Derived tables live in the Hugging Face dataset
**[conradry/biopharma-hackathon](https://huggingface.co/datasets/conradry/biopharma-hackathon)**:
the five GenomeScreen tables, plus the Parkinson's subgraph below. `data/` is gitignored, so
that repo is the durable copy.

## GenomeScreen database

`biopharma_hackathon.genomescreen` turns a [GenomeScreen](https://drug-the-whole-genome.yanyanlan.com)
release — the DrugCLIP genome-wide virtual screen — into a single DuckDB file.

```bash
uv run genomescreen-ingest /home/ubuntu/datasets/GenomeScreenDB genomescreen.duckdb
uv run genomescreen-ingest <release> <db> --limit 200   # quick smoke build
```

The 27 GB of `.pdbgz` receptor structures are left on disk; `pocket_structures` records
their paths alongside each pocket's docking grid centre.

### How the source data is keyed

- **Targets are UniProt accessions**, wrapped in AlphaFold DB ids: a directory named
  `AF-Q12879-F1-model_v4_0_pocket3` is UniProt `Q12879`, AF2 model fragment `0`, pocket `3`.
- Pockets come in two kinds. `pocketN` pockets were found by apo pocket detection;
  bare-integer pockets were transferred from an aligned holo structure, and those keep the
  source PDB id in their filenames (`..._0_5vm0_complex_refined.pdbgz` → `5vm0`).
- **Molecules are vendor catalogue ids plus SMILES — there are no PubChem CIDs.** `oid`
  is a ZINC id (`ZINC000066055208`), an Enamine REAL id (`Z1333761449_1_T2`) or an Enamine
  PV id (`PV-001914042032_1_T1`). The Enamine `_<protomer>_T<tautomer>` suffix is split out
  into `protomer_idx` / `tautomer_idx`, and `catalog_id` strips it.
- The `Name` column in `leader.csv` is *not* a molecule id — it is a row index into the
  screened library and differs between targets for the same compound. It is preserved as
  `hits.source_index`; use `catalog_id` to identify a compound.
- The same `oid` occasionally carries two SMILES (protonation variants). The `molecules`
  grain is therefore `(oid, smiles)`; group by `catalog_id` to collapse them.

### Schema

| table | grain |
| --- | --- |
| `targets` | UniProt accession |
| `pockets` | screen-result directory (target × fragment × pocket) |
| `pocket_structures` | refined receptor conformation + docking grid |
| `molecules` | distinct `(oid, smiles)` |
| `hits` | `(pocket, molecule)` DrugCLIP hit |

Two views join them up: `hit_details` (every hit with target and molecule attached) and
`target_molecule_best` (best score per target × compound, collapsing pockets).

### Querying

```python
from biopharma_hackathon.genomescreen import connect

conn = connect("genomescreen.duckdb")

# Best hits for a target
conn.sql("""
    SELECT oid, source, smiles, score, pocket_key
    FROM hit_details WHERE uniprot_acc = 'Q12879'
    ORDER BY score DESC LIMIT 10
""").show()

# Repurposing view: which targets does one compound hit?
conn.sql("""
    SELECT uniprot_acc, best_score, n_pockets
    FROM target_molecule_best WHERE catalog_id = 'ZINC000066055208'
    ORDER BY best_score DESC
""").show()

# Where to dock: grid centre for a pocket's best template structure
conn.sql("""
    SELECT template_pdb_id, structure_path, grid_center_x, grid_center_y, grid_center_z
    FROM pocket_structures WHERE pocket_key = 'AF-P12345-F1-model_v4_0_0'
""").show()
```

## Parkinson's disease subgraph

`scripts/extract_pd_subgraph.py` pulls a pathway-centric neighbourhood of Parkinson's disease
out of a [PrimeKG](https://github.com/mims-harvard/PrimeKG) edge list and writes two parquet
files. Published as the `pd_subgraph` and `pd_subgraph_enrichment` configs of the
[dataset repo](https://huggingface.co/datasets/conradry/biopharma-hackathon).

```bash
python scripts/extract_pd_subgraph.py /home/ubuntu/kg_t1.csv data/pd_subgraph.parquet
```

**PrimeKG has no `disease_pathway` relation**, so the disease reaches pathways through its
proteins: PD (`x_index` 80275) → 104 seed proteins (`disease_protein`) → 283 pathways
(`pathway_protein`) → 4,916 pathway proteins, keeping the `pathway_pathway` and
`protein_protein` edges among them. A `layer` column tags which step each edge came from.
162,384 rows, but PrimeKG stores undirected edges in both orientations — that is 81,192
unique edges mirrored, so filter `x_index < y_index` for a single copy.

Bridging through proteins drags in generic pathways alongside PD-relevant ones, so each
pathway carries a hypergeometric enrichment score against the 10,849 pathway-annotated
proteins in the graph (BH-corrected; 21 pathways reach q < 0.05, topped by catecholamine
biosynthesis, dopamine receptors and ROS detoxification). The q-value is denormalized onto
the edge table as `x_pathway_qvalue` / `y_pathway_qvalue`, populated on whichever endpoint is
a pathway. Watch out for high fold-enrichment on size-2 pathways, which are annotation
artifacts rather than real Reactome entries — add `pathway_size >= 10` for robust hits.

Two streaming passes over the source CSV, no scipy dependency: the hypergeometric survival
function and BH correction are implemented in log space in the script itself (validated
against `scipy.stats` to ~1e-13).

`scripts/build_toxin_targets.py` adds an environmental-toxin overlay on the same pathways —
ten PD-associated toxins mapped to their literature gene targets, specificity-weighted and
rolled up per pathway (`data/pd_toxin_target.parquet`, `data/pd_toxin_pathway.parquet`).
**9 of the 21 enriched pathways contain a toxin target**, topped by ROS detoxification and
mitochondrial biogenesis. Mechanism lives in the EC numbers, not the gene symbols: MPTP's
causal target (monoamine oxidase) and rotenone's (complex I) appear only as `EC:1.4.3.4` and
`EC:1.6.99.3`, so the script carries a hand-curated EC→gene map and asserts both anchors in
`tests/test_toxin_targets.py`.

## Repurposing knowledge base

`biopharma_hackathon.pdkb` joins every layer above into one DuckDB file for an agent to query,
and can unfold any candidate back into the evidence that produced it.

```bash
uv run pdkb-build                          # -> pd_kb.duckdb (~11 MB)
uv run pdkb-rank --limit 20                # the prioritised candidate list
uv run pdkb-rank --verdict mechanistic_risk
uv run pdkb-graph Rasagiline --format mermaid   # why this one might work
uv run pdkb-graph Capivasertib --format dot --out risk.dot
```

| table | rows | what it is |
| --- | ---: | --- |
| `pd_tree`, `pd_edges`, `pathway_enrichment` | 54,557 / 162,384 / 283 | PrimeKG PD neighbourhood |
| `toxin_target`, `toxin_pathway` | 79 / 76 | the toxin overlay |
| `target_drug`, `drug_summary`, `trials` | 124 / 71 / 9,802 | ChEMBL approved drugs + ClinicalTrials.gov |
| `mechanism_direction` | 30 | curated: which way to push each target |
| `screen_hits` | 6,018 | DrugCLIP hits for the 17 screened PD targets |
| `curation_conflicts` | 12 | where the two toxin curations disagree |

Views: `drug_candidates` (scored, one row per drug), `drug_target_evidence` (per drug × target),
`undrugged_targets`, `pathway_evidence`.

### Direction is the gate

PrimeKG records *that* a drug hits a protein and ChEMBL records *how*, but neither knows a toxin
and a drug can push one target opposite ways. `pdkb/mechanism.py` encodes, per target, which
action opposes the toxin insult and which mimics it. Without it the ranking promotes drugs that
reproduce the disease mechanism — it catches **venetoclax** (BCL2 inhibitor, anti-apoptotic
target), **capivasertib** (AKT inhibitor, a pathway rotenone already suppresses), **metyrosine**
(TH inhibitor, depletes dopamine) and **amphetamine** (DAT releasing agent, the route MPP+ takes
into neurons). All four score as `mechanistic_risk`, not as candidates. Targets where the
argument genuinely runs both ways — the cholinesterases — are marked `ambiguous` and scored
neutral rather than resolved by fiat.

### Ranking

`priority_score` is a weighted sum of six components, all kept as columns so it can be taken
apart: `direction` (0.30), `pathway` (0.18), `toxin` (0.18), `cns` (0.17), `clinical` (0.10),
`biomarker` (0.07). Weights live in the `score_weights` table — `UPDATE` it and re-select
`drug_candidates` to get a different ranking with no rebuild.

Validation: the three approved PD MAO-B inhibitors (rasagiline, safinamide, selegiline) come out
at ranks 2–4 without anything telling the pipeline they are PD drugs.

The `cns` component exists because trial volume alone floats systemic blockbusters to the top —
adalimumab and infliximab have 469 and 359 trials and **zero** in any neurological indication,
yet PD is a CNS disease and IgG reaches CSF at well under 1% of serum. Antibodies are discounted,
not zeroed: a peripheral immune mechanism is a live hypothesis.

### Evidence graphs

`pdkb-graph <drug>` exports the reasoning as Mermaid, Graphviz DOT or JSON:

```
toxin --implicates--> target <--acts on-- drug
                        └--member of--> pathway --enriched in--> Parkinson's
```

Drug edges are coloured by direction — green opposes the insult, red mimics it, amber is
two-sided — so a scientist sees the arguable claim rather than a number. The JSON carries the
curated rationale for each direction call.

### Caveats

- **Scores are relative to this cohort of 96 drugs**, not absolute; components are normalised
  within the set. 25 of the 96 have no trials in the source pull, so their clinical, CNS and
  biomarker components are zero by absence of data rather than by evidence — none reach the
  top 10, but check `n_trials > 0` before reading a low score as a judgement.
- **Trial counts span all indications.** Most trials are the drug's on-label use.
- **Target→drug is an annotated mechanism, not evidence of PD efficacy.** Hypothesis generation.
- **The two toxin curations disagree and the disagreement is kept.** The upstream list discards
  `nadh dehydrogenase [EC:1.6.99.3]` as a cross-species artifact, but the human NDUF subunits are
  real and rotenone's inhibition of complex I is the canonical PD mechanism. See
  `curation_conflicts`.
- **16 of 30 targets have no approved drug at all** — including SNCA, SOD2 and PGC-1β. Repurposing
  has nothing to offer them; `undrugged_targets` joins them to the virtual screen, which covers
  only 7 of the 16.
