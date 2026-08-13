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
