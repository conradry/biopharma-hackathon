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

## Data

Parkinson's drug-repurposing datasets: target genes → approved drugs → clinical
trials → outcome/biomarker measures.

```
data/
├── inputs/
│   └── pd_toxin_human_data_integrated.csv   # env-toxin → human-data → gene targets (pipeline input)
├── toxin_target/                            # from environmental-toxin gene targets (30 targets, 71 drugs)
│   ├── PD_toxin_target_Trials.csv               # main: one row per (drug, trial) + outcome measures
│   ├── PD_toxin_target_Drug_Summary.csv         # per-drug trial/biomarker counts
│   ├── PD_toxin_target_Target_Drug_map.csv      # target → approved drugs
│   ├── PD_toxin_target_Toxin_Target_map.csv     # toxin → druggable targets (+ excluded artifacts)
│   ├── PD_toxin_target_Biomarker_measures.csv   # granular biomarker outcomes (by target)
│   ├── PD_toxin_target_Biomarker_by_target.csv  # per-target biomarker rollup
│   └── PD_toxin_target_README.csv               # methods & caveats
└── pd_seed_target/                          # from PD pathway-enrichment seed genes (pd_tree.parquet; 119 drugs)
    ├── PD_seed_target_Trials.csv                # main: one row per (drug, trial) + outcome measures
    ├── PD_seed_target_Drug_Summary.csv          # per-drug trial/biomarker counts
    ├── PD_seed_target_Drug_map.csv              # drug → target
    ├── PD_seed_target_Biomarker_measures.csv    # granular biomarker outcomes (by target)
    ├── PD_seed_target_Biomarker_by_target.csv   # per-target biomarker rollup
    └── PD_seed_target_Biomarker_list.csv        # distinct biomarkers listed per target
```

All-indication scope: most trials are the drug's on-label use — filter the `Conditions`
column for Parkinson's/neuro. Target→drug is a mechanistic link (hypothesis-generating),
not evidence of repurposing efficacy.

## Skills

```
skill/
├── toxin-target-human-data/     # upstream: disease → env toxins → human data → gene targets
└── target-trial-repurposing/    # downstream: gene targets → approved drugs → trials + biomarkers
```
