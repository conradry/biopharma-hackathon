# Presentation figures

Eight PNGs at 200 dpi, one idea per figure, sized for a 16:9 slide. Regenerate with:

```bash
uv run --with matplotlib python scripts/make_presentation_figures.py             # light
uv run --with matplotlib python scripts/make_presentation_figures.py --theme both  # + figures/dark/
```

Every number is read from `pd_kb.duckdb` at render time, so re-weighting the score
(`UPDATE score_weights ...`) and re-running redraws the figures rather than leaving them stale.
The PNGs are tracked here so they can be shared and viewed directly on GitHub; they are
also always one command away if you'd rather regenerate them.

## What each one says

| file | the slide's point |
| --- | --- |
| `fig1_pipeline.png` | 10 toxins → 30 targets → 96 drugs → 44 direction-clean → 10. Every stage is a filter with a reason. |
| `fig2_toxin_convergence.png` | Chemically unrelated toxins hit the same proteins. CAT, SOD1/2, CYP1B1 and TNF are the convergence points; 16 of 30 targets have no approved drug. |
| `fig3_ranking.png` | The shortlist with its score decomposed. The bars say *why* each drug ranks, not just how high. |
| `fig4_direction_gate.png` | All 96 drugs by direction verdict. Nine push their target the way the toxin does — anti-recommendations, named. |
| `fig5_validation.png` | The positive control. Rasagiline, safinamide and methylphenidate surface at ranks 2, 4 and 6 with nothing telling the pipeline they treat PD. |
| `fig6_pathways.png` | Pathway enrichment, with convergent pathways (PD genetics *and* toxin) marked. |
| `fig7_mechanism_network.png` | Toxin → target → drug. The ten names are three mechanisms: 5 anti-TNF, 4 MAO, 1 DAT. |
| `fig8_biomarkers.png` | 76% of annotated trial measures are pharmacokinetics. Excluding them is what makes the biomarker component mean "testable". |

## Suggested deck order

`fig1` (what we searched) → `fig2` (the toxin signal) → `fig6` (pathway corroboration) →
`fig4` (the gate) → `fig5` (validation) → `fig3` (the answer) → `fig7` (correlated risk) →
`fig8` (how you'd test it).

If you only get three slides: **fig2, fig5, fig3**.

## Caveats to keep on the slide

- Scores are normalised across the 96-drug cohort. 0.77 means "top of this list", not a
  probability of success. Quote the components, not the total.
- 25 of the 96 drugs have no trials in the source pull, so their clinical, CNS and biomarker
  components are zero from missing data, not from a judgement.
- Nothing here is efficacy data — it is mechanistic annotation plus trial metadata.

See `data/README.md` for how the ranking itself is built and what the columns mean.

## Colour

Both themes are validated categorical palettes (six slots, adjacent-pair CVD ΔE ≥ 8, normal-vision
ΔE ≥ 15) against their own surface; the dark set is re-stepped for the dark surface rather than
flipped. Sequential encodings use a single blue ramp. Values are directly labelled everywhere,
which is what discharges the sub-3:1 contrast on the lighter categorical slots in light mode.
