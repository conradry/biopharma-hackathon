---
name: pd-drug-repurposing
description: Answer questions about repurposing approved drugs for Parkinson's disease from the pd_kb.duckdb evidence base — which candidates look most promising and why, what the toxin/pathway/trial evidence says, which drugs are mechanistically wrong for PD, which targets have no drug at all, and exportable evidence graphs that explain a hypothesis to a scientist. Use this skill whenever the user asks about promising or candidate drugs for Parkinson's, PD drug repurposing, which compounds to prioritize or avoid, why a particular drug might work for PD, PD targets or pathways, environmental-toxin evidence for PD, or wants a repurposing rationale written up — even if they never mention the database, DuckDB, or the file by name.
---

# Parkinson's drug repurposing

You have an assembled evidence base that ranks approved drugs as Parkinson's repurposing
candidates and can unfold any of them back into the evidence that produced the ranking.

Your job is almost never "run the query and paste the list." The ranking is a starting point
that needs interpretation, and the interpretation is where the value is. A scientist reading
your answer should come away knowing what to try, what to avoid, and what would change their
mind.

## Getting to the data

The database lives at `pd_kb.duckdb` in the repo root. It is gitignored, so it may not exist:

```bash
ls pd_kb.duckdb || uv run pdkb-build     # ~30s, needs data/ populated
```

Query it read-only. DuckDB's CLI or Python both work:

```bash
uv run python -c "
import duckdb; c = duckdb.connect('pd_kb.duckdb', read_only=True)
print(c.execute('SELECT drug, verdict, priority_score FROM drug_candidates LIMIT 5').fetchall())
"
```

Two CLIs wrap the common cases: `uv run pdkb-rank [--limit N] [--verdict V] [--novel-only]`
and `uv run pdkb-graph <drug> --format mermaid|dot|json`.

## The one thing to get right

**A toxin and a drug can hit the same target in opposite directions, and only one of those
directions is a therapy.** The `direction` column carries this, and it changes the conclusion
rather than nudging a score:

- `protective` — the drug opposes the toxic insult. This is the hypothesis.
- `risk` — the drug pushes the target the same way the toxin does. Capivasertib inhibits AKT1,
  a survival pathway rotenone already suppresses. Metyrosine inhibits tyrosine hydroxylase in a
  disease defined by failing dopamine synthesis. These are not weak candidates; they are the
  wrong hypothesis, and proposing one would be a serious error.
- `ambiguous` — the argument genuinely runs both ways and the curation refuses to resolve it.
  The cholinesterases are the live example: chlorpyrifos causes harm by inhibiting AChE, yet
  donepezil and rivastigmine are inhibitors used in PD dementia. Present both sides.
- `unknown` — no defensible direction. Never report this as support.

So `verdict = 'mechanistic_risk'` drugs are excluded from proposals, and worth mentioning
explicitly when the user asks what to avoid or when one would otherwise have scored well.

## Separating proposals from validation

This trips up every naive pass. The top of the raw ranking is dominated by drugs that are
*already* PD therapies — rasagiline, safinamide, selegiline, rivastigmine, donepezil,
methylphenidate. Their presence is the pipeline validating itself: nothing told it these were
PD drugs and it surfaced them anyway. That is worth reporting as a confidence check.

It is **not** an answer to "what should we try that's new." For proposals, exclude them:

```sql
SELECT drug, targets, modality, priority_score, direction_component,
       pathway_component, toxin_component, cns_component, n_trials, n_neuro_trials
FROM drug_candidates
WHERE verdict = 'candidate' AND n_pd_trials = 0
ORDER BY priority_score DESC;
```

Lead with the validation in a sentence, then spend the answer on the novel candidates.

**Diversify across mechanisms, not just rows.** The top of that list is currently four anti-TNF
antibodies — etanercept, adalimumab, infliximab, certolizumab. Handing back three of them looks
like three suggestions but is one hypothesis with a shared failure mode: if TNF inhibition
doesn't reach the brain, all four die together. When a user asks for a handful of candidates,
give them distinct targets or mechanisms, and say plainly when the list collapses onto one
idea. `GROUP BY targets` or scan the `targets` column before choosing.

## Reading a score honestly

`priority_score` is a weighted sum of six components, all present as columns:
`direction` (0.30), `pathway` (0.18), `toxin` (0.18), `cns` (0.17), `clinical` (0.10),
`biomarker` (0.07). Weights live in the `score_weights` table — an agent can `UPDATE` it and
re-select `drug_candidates` for a different ranking with no rebuild, which is the right move
when a user says something like "I care more about mechanism than trial history."

Three things to keep straight when you describe a score:

- **It is relative to this cohort of 71 drugs**, normalised within the set. 0.76 is "top of
  this list," not a probability of success. Say it that way.
- **Report the components, not just the total.** Two drugs at 0.65 can be there for completely
  different reasons, and which reason it is determines whether the scientist believes it.
- **A high score is a hypothesis worth testing, not evidence of efficacy.** Everything here is
  mechanistic annotation plus trial metadata. No efficacy data for PD is in this database.

## The main tables

Views first — they carry the joins you want.

| view | grain | use it for |
| --- | --- | --- |
| `drug_candidates` | drug (71) | the ranking, components, verdict, modality |
| `drug_target_evidence` | drug × target | direction, rationale, per-target evidence |
| `undrugged_targets` | target (16) | PD targets no approved drug engages |
| `pathway_evidence` | pathway (283) | where the two evidence routes converge |

Underlying tables: `pd_tree` (disease→pathway→protein→drug), `pathway_enrichment`,
`toxin_target` / `toxin_pathway` (the toxin overlay), `target_drug` / `drug_summary` /
`trials` (ChEMBL + ClinicalTrials.gov, 9,802 trials), `mechanism_direction` (the curation,
with a `rationale` column worth quoting), `screen_hits` (DrugCLIP virtual screen),
`curation_conflicts`.

`verdict` values: `candidate`, `established_pd_therapy`, `mechanistic_risk`,
`direction_ambiguous`, `direction_unknown`.

## Common question shapes

**"What should we repurpose for PD?"** — Validation sentence, then the novel-candidate query
above. For each of the top few, give the target, the direction rationale, the toxin and pathway
evidence, and the specific reason it might fail. Three well-explained candidates beat ten rows.

**"Why might drug X work?"** — `pdkb-graph X --format mermaid` for the picture, plus:

```sql
SELECT gene, action_type, direction, direction_confidence, direction_rationale,
       toxins, best_pathway_qvalue, n_enriched_pathways
FROM drug_target_evidence WHERE lower(drug) = lower('X');
```

The `direction_rationale` is written to be quotable. Include the trials that already exist:
`SELECT nct_id, phase, status, conditions FROM trials WHERE lower(drug)=lower('X') AND is_neuro`.

**"What should we avoid?"** — `WHERE verdict = 'mechanistic_risk'`, with the reason from
`mechanism_direction.rationale` for each.

**"Which targets need new chemistry?"** — `undrugged_targets`. 16 of 30 PD-implicated targets
have no approved drug by direct mechanism, including SNCA, SOD2 and PGC-1β. Repurposing
structurally cannot serve these, which is often the most useful thing to tell a user chasing
disease modification. `n_screen_hits` shows whether the virtual screen offers a starting point
(it covers only 7 of the 16, and not the interesting ones).

**"Which pathways?"** — `pathway_evidence WHERE is_convergent` gives the 11 pathways that are
both enriched for PD proteins and hit by a toxin. ROS detoxification and mitochondrial
biogenesis top it, which is the canonical PD mechanism arriving from two independent directions.

## Evidence graphs

When the user wants to explain a hypothesis to a scientist, export the graph:

```bash
uv run pdkb-graph Etanercept --format mermaid        # inline in markdown
uv run pdkb-graph Etanercept --format dot --out e.dot   # publication figure
uv run pdkb-graph Etanercept --format json           # includes direction rationale
```

It renders `toxin → target ← drug`, `target → pathway → Parkinson's`, plus neuro trials, with
drug edges coloured by direction. Paste the Mermaid straight into your answer — the coloured
direction edge is the claim the scientist will want to argue with, which is the point.

## Caveats to surface

Not all of these in every answer — pick the ones that bear on the question. But do not let a
user walk away believing the coverage is broader than it is.

- **This is a narrow slice, not a comprehensive PD screen.** 71 drugs, reached from 30 targets,
  reached from 10 environmental toxins. Plenty of good PD repurposing hypotheses are simply not
  in here. Say so when a user asks "what are the best candidates" as though the list were
  exhaustive.
- **Trial counts span every indication.** 469 adalimumab trials are mostly rheumatology. Trial
  volume is a safety-data proxy, not PD evidence.
- **Antibodies are discounted for a CNS disease and the discount is a judgement call.** The
  anti-TNF drugs have strong target-level evidence and near-zero brain penetration; a peripheral
  immune mechanism is a live hypothesis, not a settled one. Flag `modality = 'biologic'` when it
  appears near the top.
- **Target→drug is an annotated mechanism, not evidence the drug was ever tried in PD.**
- **The toxin evidence is literature co-mention, not affinity**, and the source lists only the
  top 8 targets per toxin.
- **Two curations disagree about complex I.** `curation_conflicts` records it: the upstream
  target list discarded `nadh dehydrogenase` as a cross-species artifact, but the human NDUF
  subunits are real and rotenone's complex I inhibition is the canonical PD mechanism. It has no
  approved drug, so it contributes nothing to the ranking — worth raising when a user asks about
  mitochondrial mechanisms or why rotenone's main target is missing.

## Writing the answer

Lead with the finding, not the method. A scientist wants the candidate and the reason before
they want your query.

For each candidate you propose, cover: **target and action**, **why that direction is the
protective one**, **the toxin and pathway evidence**, **clinical readiness**, and **the specific
reason it might fail**. That last one is not hedging — the failure mode is usually the most
informative thing you can say, and omitting it makes the answer less useful, not more
confident.

Quote real numbers (q-values, trial counts, toxin counts) rather than adjectives. Offer the
evidence graph when the user is going to show this to someone else.
