# Skills

Two portable [Claude Code / Paperclip](https://claude.com/claude-code) skills that produce the
Parkinson's datasets in `data/`. Each is a self-contained folder (`SKILL.md` + `scripts/`);
copy one into `.claude/skills/` (or a Paperclip skills dir) to run it. Both are disease-agnostic
— swap the disease/target inputs.

| skill | what it does | feeds |
| --- | --- | --- |
| [`toxin-target-human-data/`](toxin-target-human-data/) | disease → associated environmental toxins (PubMed) → human-relevant safety data + literature gene targets (PubChem) → one integrated CSV | `data/pd_toxin_human_data_integrated.csv` |
| [`target-trial-repurposing/`](target-trial-repurposing/) | gene targets → UniProt → approved drugs (ChEMBL) → all-indication clinical trials + outcome/biomarker measures (ClinicalTrials.gov), organised by target | `data/toxin_target/*`, `data/pd_seed_target/*` |

They chain: the first turns a disease into a toxin→target table, the second turns targets into a
drug→trial→biomarker evidence base.

## target-trial-repurposing

Runs on the Paperclip CLI (ChEMBL + UniProt + ClinicalTrials.gov/AACT). Four steps, each a script:

```bash
python3 scripts/derive_drugs.py         --genes SNCA,MAOB,DRD2 --mode direct --out drugs.tsv
python3 scripts/extract_target_trials.py --drugs drugs.tsv --out trials.tsv
python3 scripts/extract_biomarkers.py    --drugs drugs.tsv --out biomarkers.tsv
python3 scripts/assemble_xlsx.py         --trials trials.tsv --out db.xlsx
python3 scripts/assemble_biomarkers.py   --in biomarkers.tsv --out-dir data --prefix "myproj_"
```

The extraction scripts work around two Paperclip quirks that make a naive `SELECT` lossy — a
~60-char per-cell display cap (long outcome text is chunked and reassembled) and a ~48KB stdout
cap (queries are paged at 60 rows). See each script's docstring and `target-trial-repurposing/SKILL.md`.

## toxin-target-human-data

Python + `requests` only (PubMed E-utilities, PubChem PUG REST/View). See its `SKILL.md` and
`reference/parkinsons_example.md` for the worked Parkinson's example.
