# PD-seed pathway repurposing dataset

A target→approved-drug→clinical-trial evidence base for Parkinson's, built from the
**PD-seed genes** of the pathway-enrichment tree (`pd_tree.parquet`) rather than the
environmental-toxin targets. It is a parallel, broader cut to the toxin pipeline that feeds
`pdkb`: same method, a different (larger, more dopaminergic) target set.

## How it was built

```
pd_tree.parquet  →  90 PD-seed genes (protein_is_pd_seed)
                 →  UniProt (Homo sapiens)
                 →  119 approved drugs   (ChEMBL drugs_by_accession, max_phase = 4, direct mechanism)
                 →  8,940 trials         (ClinicalTrials.gov / AACT, all indications)
                 →  outcome + biomarker measures, organised by target
```

Approved = ChEMBL `max_phase = 4`. Drugs are matched to trials on the intervention name
(`interventions.name ILIKE '%drug%'`), so all salt forms of a drug are captured. Built with
the `target-trial-repurposing` skill (see `skill/`).

The target set is core PD/CNS pharmacology, which is why it separates cleanly from the toxin
cut: **DRD2** (dopamine agonists + antipsychotics, 59 drugs), **SLC6A3/DAT**, **MAO-A/B**,
**SLC18A2/VMAT2** (tetrabenazine class), **DRD1**, plus TNF, IGF1R, INSR, BCHE, DDC (carbidopa),
TH (metyrosine). 571 of the 8,940 trials are in a PD/neuro indication — ~6× the toxin cut's rate.

## Files

| file | grain | rows |
| --- | --- | ---: |
| `PD_seed_target_Drug_map.csv` | drug → target(s) | 119 |
| `PD_seed_target_Drug_Summary.csv` | one row per drug (with trials): `N_trials`, `N_biomarker_trials`, `N_PDneuro_trials` | 91 |
| `PD_seed_target_Trials.csv` | one row per (drug, trial): conditions, phase, status, primary/secondary outcome measures, biomarker flag | 8,940 |
| `PD_seed_target_Biomarker_measures.csv` | one row per (target, trial, biomarker outcome): full measure text, category, time frame, PD/neuro flag | 4,149 |
| `PD_seed_target_Biomarker_by_target.csv` | per-target rollup: measure counts by category + distinct measures | 17 |
| `PD_seed_target_Biomarker_list.csv` | distinct (target, biomarker measure) with trial counts + the drugs measuring it | 3,600 |

Only 91 of the 119 drugs have any trial; only 17 of the targets carry a biomarker-bearing
trial (many PD-seed genes — SNCA, LRRK2, PRKN, PINK1, VPS35 — have no approved direct-mechanism
drug at all).

`Biomarker_category` ∈ `PK · imaging · fluid_molecular · inflammation · genomic · other_biomarker`,
assigned by priority-ordered keyword match on the measure text. The disease-modification-relevant
readouts are imaging/fluid_molecular/inflammation — e.g. carbidopa's **[18F]-FDOPA PET**,
rotigotine's **plasma neurofilament**, bromocriptine's **CSF Aβ / tau-PET**. PK dominates the raw
counts and only tells you the drug reached plasma; filter it out and require `PD_neuro = Yes` for
the interesting subset.

## Caveats

- **All indications.** Most trials are the drug's on-label use (antipsychotic, oncology,
  endocrine), not Parkinson's — filter the `Conditions` column.
- **Target→drug is an annotated mechanism, not evidence of PD efficacy.** Hypothesis generation.
  It also says nothing about *direction*: many of these drugs (DRD2 antagonists, VMAT2 inhibitors,
  metyrosine) push dopaminergic targets the *wrong* way for PD. The `pdkb` mechanism gate handles
  that for the toxin cut; this dataset does not encode it.
- **Outcome text is capped (~220 chars)** in `Trials.csv`; `N_primary`/`N_secondary` give the true
  counts and the `NCT_ID`/`URL` link the full record. The biomarker files carry the full measure text.
