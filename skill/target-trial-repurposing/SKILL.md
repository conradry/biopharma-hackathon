---
name: target-trial-repurposing
description: Drug-repurposing pipeline over Paperclip. Given a set of gene targets (typed, or a column in a CSV such as an enrichment/toxin table), derive the approved drugs that act on those targets (ChEMBL), then pull ALL of those drugs' clinical trials across every indication from ClinicalTrials.gov with per-trial outcome measures and a biomarker flag, and assemble a styled multi-sheet Excel database. Use when a user wants "clinical data / outcome measures / biomarker measures for approved drugs that hit these targets", target->drug->trial mapping, or a repurposing evidence base.
---

# Target → Approved-Drug → Clinical-Trial repurposing pipeline

Turns a list of **gene targets** into a structured **clinical-trial outcome database** for
the approved drugs that hit those targets. Runs on the `paperclip` CLI (ChEMBL + UniProt +
ClinicalTrials.gov/AACT). Load the Paperclip skill first (`paperclip skill`) if unfamiliar.

## When to use
- "Find clinical trials / outcome measures / biomarkers for approved drugs targeting <genes>."
- A user has a target list (or a CSV column of targets, e.g. a pathway/toxin enrichment table)
  and wants the matching approved drugs' trial evidence for repurposing.

## Inputs you need from the user (ask if unclear)
1. **Targets** — gene symbols, or a CSV + column name. Extract symbols from bracketed tokens
   like `... [SNCA] (868)`; **exclude non-human ortholog artifacts** common in KEGG/BlastKOALA
   output (viral e.g. influenza PB2, bacterial e.g. mercuric reductase, plant/fungal e.g.
   laccase/peroxidase/photosystem II). Map ambiguous enzyme names to human genes
   (monoamine oxidase→MAOA/MAOB, tyrosine hydroxylase→TH, acetylcholinesterase→ACHE, etc.).
2. **Derivation mode** — `direct` (annotated mechanism, well-drugged targets) or `bioactivity`
   (approved compounds with measured potency; wider repurposing net, noisier).
3. **Trial scope** — all indications (default here) vs. filter to a disease. Warn the user that
   "all indications" for well-drugged targets is large (well-drugged drugs have 100s–1000s of
   trials each, mostly on-label). Give counts before committing.
4. **Output** — xlsx (default), CSV also emitted.

## Prerequisites
- `paperclip config` authenticated. `python3 -m pip install openpyxl`.
- Before ANY protein SQL, read `paperclip skills show proteins` (schema for UniProt/ChEMBL views).

## Workflow

**Step 1 — targets → UniProt → approved drugs → `drugs.tsv`**
```
python3 scripts/derive_drugs.py --genes SNCA,MAOB,ACHE,KDR --mode direct --out drugs.tsv
```
`drugs.tsv` columns: `search_term<TAB>display<TAB>targets`. **Review it** and delete any
over-generic `search_term` (e.g. `SYNTHETIC`, `ESTROGENS`) that would over-match trial
interventions. Salt forms are already collapsed to the base ingredient.

Sanity-check trial volume before extracting (per drug):
```
paperclip sql -s trials "SELECT COUNT(DISTINCT i.nct_id) FROM ctgov.interventions i WHERE i.name ILIKE '%DONEPEZIL%'"
```

**Step 2 — drugs → trials + outcomes → `trials_rows.tsv`**
```
python3 scripts/extract_target_trials.py --drugs drugs.tsv --out trials_rows.tsv
```
Run it in the background for large sets; it prints per-drug progress. One row per (drug, trial):
targets, drug, NCT, phase, status, start_date, enrollment, conditions, n_primary, n_secondary,
primary/secondary outcome-measure text, biomarker flag.

**Step 3 — assemble the Excel database**
```
python3 scripts/assemble_xlsx.py --trials trials_rows.tsv --out database.xlsx --title "..."
```
Sheets: README (methods+caveats), Drug_Summary (per-drug trial/biomarker counts), Trials (main).
To merge batches, pass multiple TSVs: `--trials batchA.tsv batchB.tsv` (deduped on drug+NCT).

**Step 4 (optional) — biomarker measures organized BY TARGET**
The Trials table gives a coarse per-trial biomarker Yes/No. For an itemized, target-centric
view (what biomarker endpoints exist across each target's drug trials):
```
python3 scripts/extract_biomarkers.py --drugs drugs.tsv --out biomarker_rows.tsv
python3 scripts/assemble_biomarkers.py --in biomarker_rows.tsv --out-dir data --prefix "myproj_"
```
Produces:
- `<prefix>Biomarker_measures.csv` — granular, one row per (target × trial × biomarker outcome
  measure): full measure text, `Biomarker_category` (PK/imaging/fluid_molecular/inflammation/
  genomic/other_biomarker), time frame, condition, PD/neuro flag.
- `<prefix>Biomarker_by_target.csv` — per-target rollup: #drugs, #trials, #measures, category
  counts, and the distinct measures aggregated under each target.

`extract_biomarkers.py` filters `ctgov.design_outcomes` to biomarker-matching measures (edit `BT`
to tune the net) and explodes each drug's outcomes across all targets it hits. Categorization is
keyword-based and priority-ordered so a specific analyte beats generic PK; keep keywords specific
(e.g. `" spect "` not `"spect"`, which would hit "inspection"; `"tau protein"` not `" tau "`,
which hits the PK "tau interval"). Category is a convenience — the full `Measure` text is authoritative.
Most all-indication biomarkers are PK; the repurposing-interesting rows are imaging/fluid_molecular/
inflammation — filter on `Biomarker_category` and `PD_neuro`.

## Why the extraction script is non-trivial (do not "simplify" it away)
Paperclip's trials SQL is delivered through an ASCII renderer with two hard limits the script
works around, so a naive `SELECT ... measure` loses data:
- **~60-char per-cell cap** → long outcome text is chunked into 55-char `SUBSTRING` columns,
  sentinel-padded (`¤`) so renderer spaces can't corrupt boundaries, and reassembled in Python.
- **~48KB stdout cap** → page at **60 wide rows/query** (`LIMIT/OFFSET`); larger pages silently
  truncate. Narrow queries can go to ~200.
- Free-text measures contain newlines/pipes → `translate()`'d out in SQL so every row is one
  pipe-safe line.
- The biomarker flag is computed in SQL with `bool_or(measure ILIKE ANY(...))` over ALL of a
  trial's outcomes (untruncated) — more complete than matching the capped text. Edit
  `BIOMARKER_TERMS` in `extract_target_trials.py` to tune the net; keep patterns specific
  (avoid `'% spect%'`, which hits "respect/aspect").

## Key Paperclip facts encoded here
- Approved = ChEMBL `max_phase=4`. `drugs_by_accession` = direct mechanism only;
  `bioactivities_by_accession` = any measured activity (has `pchembl_value`).
- Join key across UniProt/PDB/ChEMBL is the UniProt `accession`.
- AACT tables used: `ctgov.studies`, `ctgov.interventions`, `ctgov.conditions`,
  `ctgov.design_outcomes` (`outcome_type` ∈ primary/secondary/other).
- chembl/proteins SQL previews large results (page at LIMIT 25); trials SQL dumps rows.

## Caveats to surface in every deliverable
All-indication scope means most trials are on-label, not the target disease — filter Conditions.
Target→drug is a mechanistic link, not evidence of repurposing efficacy. Outcome text is capped;
counts + NCT link give the full picture. Biomarker flag is a screening aid to verify per trial.
