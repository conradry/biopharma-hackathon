# Parkinson's repurposing outputs

Generated artifacts: a ranked shortlist of repurposing candidates and the graphs that
explain it. Everything here is derived — regenerate with:

```bash
uv run pdkb-build                              # data/ + HF -> pd_kb.duckdb
uv run python scripts/export_ranking_graph.py  # pd_kb.duckdb -> the files below
```

The database itself (`pd_kb.duckdb`, ~11 MB) is gitignored; these small outputs are tracked.

## Files

| file | what it is |
| --- | --- |
| `pd_ranking_top10.csv` | the shortlist, one row per drug, rank 1–10 |
| `pd_top10_nodes.csv` / `pd_top10_edges.csv` | mechanism graph as tables (Gephi, Cytoscape, pandas) |
| `pd_top10.graphml` | the same graph in one file |
| `pd_top10_mechanism.mmd` | Mermaid: why these drugs are on the list (24 nodes, 55 edges) |
| `pd_top10_readouts.mmd` | Mermaid: how you'd measure whether one works (12 nodes, 22 edges) |
| `pd_pathway_convergence.mmd` | Mermaid: pathway-level evidence (41 nodes, 77 edges) |

Every graph is capped at 100 edges. The exporter reports anything it drops rather than
truncating silently — all three currently fit with room to spare.

Paste a `.mmd` file straight into any Markdown that renders Mermaid. For `.graphml`, open in
Cytoscape or Gephi; node and edge attributes (scores, q-values, direction) come through as
columns you can style on.

## Reading the ranking

One row per drug, ranked by `priority_score`. The columns that matter most:

- **`verdict`** — all ten are `candidate`. The other values in the database are
  `established_pd_therapy`, `mechanistic_risk`, `direction_ambiguous`, `direction_unknown`.
- **`direction_component`** — 1.0 means the drug opposes the toxic insult at every target it
  hits. This is the gate, not a tiebreaker; see below.
- **`modality`** — `biologic` flags a large molecule that does not meaningfully cross the
  blood-brain barrier.
- **`cns_component`** — 1.0 existing PD trials, 0.7 other neuro trials, 0.15 biologic with
  neither, 0.5 small molecule with neither.
- **`biomarker_component`** — share of trials with a *disease-relevant* readout (imaging, fluid
  molecular, inflammation, genomic). Pharmacokinetics is excluded on purpose: it is 76% of the
  annotated measures and only tells you the drug reached plasma.
- **`n_pd_trials`** is 0 for every row by construction — these are proposals, not things already
  being tried.

`priority_score` = 0.30·direction + 0.18·pathway + 0.18·toxin + 0.17·cns + 0.10·clinical +
0.07·biomarker. Weights live in the `score_weights` table; `UPDATE` it and re-select
`drug_candidates` to re-rank without rebuilding.

## What another agent should know before using this

**The score is relative, not absolute.** Components are normalised across the 96-drug cohort.
0.77 means "top of this list", not a probability of success. Report the components, not the
total — two drugs at 0.65 can be there for entirely different reasons.

**Direction changes the conclusion, not the score.** A toxin and a drug can hit one target in
opposite ways. Drugs marked `mechanistic_risk` in the database (venetoclax, capivasertib,
metyrosine, amphetamine) push targets the same way the toxin does; they are the wrong
hypothesis, not weak candidates, and are excluded from this shortlist.

**The top 10 is three mechanisms, not ten ideas.** Five anti-TNF biologics, four MAO
inhibitors, one dopamine transporter inhibitor. If TNF inhibition fails to reach the brain,
five of these die together. Treat the list as three bets with redundancy, and pick across
groups rather than down the ranking.

**Already-approved PD drugs were excluded, and that they scored highly is the validation.**
Rasagiline, safinamide and selegiline rank 2–4 in the unfiltered ranking with nothing telling
the pipeline they are PD drugs. They appear in the mechanism graph as dashed anchors.

**25 of the 96 drugs have no trials in the source pull.** Their clinical, CNS and biomarker
components are zero because data is missing, not because the drug is unpromising. None reach
the top 10, but check `n_trials > 0` before reading a low score as a judgement.

**This is a narrow slice.** 96 drugs, reached from 30 targets, reached from 10 environmental
toxins. Good PD repurposing hypotheses that don't touch those 30 targets are simply absent.
Nothing here is efficacy data — it is mechanistic annotation plus trial metadata.

## Going deeper

The shortlist is a view over `pd_kb.duckdb`, which holds the full evidence base — PrimeKG
neighbourhood, pathway enrichment, toxin overlay, ChEMBL mechanisms, 9,802 trials, 6,953
categorised biomarker measures, and the DrugCLIP virtual screen. Useful entry points:

```sql
SELECT * FROM drug_candidates ORDER BY priority_score DESC;   -- all 96
SELECT * FROM drug_target_evidence WHERE drug = 'Etanercept'; -- per-target, with rationale
SELECT * FROM undrugged_targets;      -- 16 PD targets no approved drug engages
SELECT * FROM pathway_evidence WHERE is_convergent;  -- 11 pathways, both evidence routes
SELECT * FROM drug_biomarker_readouts WHERE drug = 'Etanercept';
```

For a single drug's evidence graph: `uv run pdkb-graph <drug> --format mermaid|dot|json`.

The `pd-drug-repurposing` skill in `.claude/skills/` carries the full interpretation guide and
the query recipes for common questions.
