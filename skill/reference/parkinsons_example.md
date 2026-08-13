# Worked example — Parkinson's disease environmental toxins

Ran the full pipeline with the default PD query over 100 PubMed papers. 18 distinct toxins
surfaced; 10 resolved to a PubChem compound with human-relevant data and gene targets. The
8 unresolved entries are exposure classes / abbreviations (air pollution, particulate matter,
PM2.5, traffic-related, TCE, organochlorine, organophosphate, polychlorinated biphenyl) with
no single CID.

Result file: `../../data/pd_toxin_human_data_integrated.csv`.

## Resolved toxins and their top gene targets

| Toxin | CID | PubMed mentions | Top gene targets (symbol) |
|-------|-----|-----------------|---------------------------|
| paraquat | 15939 | 14 | SOD (EC:1.15.1.1), CAT, TNF, SOD2 |
| lead | 5352425 | 12 | PB2, CAT, SOD, CD2 |
| rotenone | 6758 | 7 | **SNCA**, **TH**, NADH dehydrogenase, CASP3 |
| MPTP | 1388 | 6 | **TH**, **SNCA**, monoamine oxidase, **SLC6A3** |
| trichloroethylene | 6575 | 5 | LON, CYP1B1, CYP2E1 |
| manganese | 23930 | 3 | SOD, SOD2, CAT |
| mercury | 23931 | 3 | CD2, CAT, SOD, glutathione peroxidase |
| chlorpyrifos | 2730 | 2 | **ACHE**, BCHE, CAT |
| tetrachloroethylene | 31373 | 1 | LON, PPARGC1B, CYP |
| dieldrin | 969491 | 1 | FER1HCH, CYP1B1, ACHE |

## Why the targets are a useful signal

The classic PD neurotoxins **MPTP** and **rotenone** both top out on **tyrosine hydroxylase
(TH)** and **alpha-synuclein (SNCA)** — the core dopaminergic / Lewy-body biology of
Parkinson's — and MPTP additionally hits the dopamine transporter **SLC6A3 (DAT)**, its known
uptake route. **Chlorpyrifos**, an organophosphate pesticide, correctly resolves to
**acetylcholinesterase (ACHE)**. These are not hand-curated; they fall straight out of the
literature co-occurrence data, which is what makes the target column a fast plausibility
filter for a repurposing / mechanism hypothesis.
