"""Which direction of target engagement is plausibly protective in Parkinson's.

PrimeKG records that a drug acts on a protein but never how, and ChEMBL's
``action_type`` records how but not whether that helps. Neither knows that a
toxin and a drug hitting the same target can push it opposite ways. Without this
table a ranking treats *any* engagement of a PD-implicated target as good, which
promotes drugs that reproduce the toxic insult -- an AKT inhibitor for a pathway
that rotenone already suppresses, a tyrosine hydroxylase inhibitor for a disease
defined by failing dopamine synthesis.

So each target carries the actions that oppose the insult, the actions that
mimic it, and how firmly that is known. ``confidence`` is the honest part:

    established  approved PD drugs or human evidence support the direction
    plausible    mechanism is clear and supported in models, untested in PD
    ambiguous    defensible arguments both ways -- do not let a score decide
    unknown      no defensible direction; scored as neutral, never as support

Nothing here is a claim of efficacy. It only says which way to push a target if
you were going to push it, and refuses to guess when that is genuinely unclear.
"""

from __future__ import annotations

from dataclasses import dataclass

# ChEMBL action_type values in the source: INHIBITOR, AGONIST, ANTAGONIST,
# MODULATOR, ACTIVATOR, DEGRADER, ANTAGONIST/DEGRADER, RELEASING AGENT.
NEGATIVE = ("INHIBITOR", "ANTAGONIST", "DEGRADER", "ANTAGONIST/DEGRADER")
POSITIVE = ("AGONIST", "ACTIVATOR")


@dataclass(frozen=True)
class Direction:
    """How a target should be pushed, and how confident that is."""

    gene: str
    protective: tuple[str, ...]
    risk: tuple[str, ...]
    confidence: str
    rationale: str


DIRECTIONS: tuple[Direction, ...] = (
    Direction(
        "MAOB",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="established",
        rationale=(
            "MAO-B converts MPTP to the neurotoxic MPP+, and its oxidation of dopamine "
            "generates hydrogen peroxide. Selegiline, rasagiline and safinamide are "
            "approved for PD; MAO-B inhibition blocks MPTP toxicity in models."
        ),
    ),
    Direction(
        "MAOA",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="plausible",
        rationale=(
            "Shares the oxidative-deamination mechanism with MAO-B, but MAO-A inhibition "
            "carries tyramine/hypertensive-crisis risk and is not the PD-selective arm. "
            "Treat MAO-A-only inhibitors as weaker than MAO-B-selective ones."
        ),
    ),
    Direction(
        "SLC6A3",
        protective=("INHIBITOR", "ANTAGONIST"),
        risk=("RELEASING AGENT", "AGONIST"),
        confidence="plausible",
        rationale=(
            "The dopamine transporter is how MPP+ enters dopaminergic neurons, so DAT "
            "blockade is protective against the toxin. Substrates and releasing agents "
            "(amphetamine) do the opposite while sharing the same target annotation -- "
            "this is the clearest case where action_type flips the interpretation."
        ),
    ),
    Direction(
        "AKT1",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale=(
            "PI3K/AKT signalling is pro-survival in dopaminergic neurons and is suppressed "
            "by rotenone and MPTP. Oncology AKT inhibitors therefore push the same "
            "direction as the toxin and are flagged as a mechanistic risk, not a candidate."
        ),
    ),
    Direction(
        "BCL2",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale=(
            "Anti-apoptotic. Toxin models kill dopaminergic neurons through the intrinsic "
            "apoptotic pathway, so BCL2 inhibitors (venetoclax class) oppose the "
            "protective direction."
        ),
    ),
    Direction(
        "CASP3",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="plausible",
        rationale=(
            "Executioner caspase downstream of the same apoptotic cascade; inhibition is "
            "protective in toxin models. No approved selective inhibitor exists."
        ),
    ),
    Direction(
        "TNF",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="plausible",
        rationale=(
            "Microglial TNF drives the neuroinflammatory component of toxin models, and "
            "epidemiology of anti-TNF therapy in inflammatory disease reports lower PD "
            "incidence. CNS exposure of biologics is the limiting question, not direction."
        ),
    ),
    Direction(
        "IL6",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="plausible",
        rationale=(
            "Same neuroinflammatory rationale as TNF, on weaker evidence. Elevated in PD "
            "serum and CSF, but causality is much less established."
        ),
    ),
    Direction(
        "CAT",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale=(
            "Catalase clears the hydrogen peroxide that paraquat, rotenone and MPTP "
            "generate. Augmenting antioxidant capacity is protective; inhibiting it "
            "reproduces the insult. No approved activator exists."
        ),
    ),
    Direction(
        "SOD1",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale="Superoxide dismutation is the first step of the same detoxification arm.",
    ),
    Direction(
        "SOD2",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale=(
            "The mitochondrial isoform, directly downstream of complex I superoxide "
            "leak -- the most mechanistically apt member of the antioxidant set."
        ),
    ),
    Direction(
        "GPX1",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale="Glutathione peroxidase clears peroxides; glutathione is depleted in PD nigra.",
    ),
    Direction(
        "PPARGC1B",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale=(
            "The PGC-1 axis drives mitochondrial biogenesis and is downregulated in PD "
            "substantia nigra; activation opposes the complex I deficit."
        ),
    ),
    Direction(
        "GABPA",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="plausible",
        rationale="NRF-2/GABP is the transcription factor arm of the same biogenesis programme.",
    ),
    Direction(
        "SNCA",
        protective=("DEGRADER", "ANTAGONIST/DEGRADER", "INHIBITOR"),
        risk=POSITIVE,
        confidence="established",
        rationale=(
            "Lowering alpha-synuclein is the central disease-modifying hypothesis in PD; "
            "SNCA multiplication causes familial disease. Direction is not in doubt, "
            "though no approved drug engages it directly."
        ),
    ),
    Direction(
        "TH",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="established",
        rationale=(
            "Tyrosine hydroxylase is rate-limiting for dopamine synthesis, and its loss "
            "defines the disease. Inhibitors (metyrosine) deplete dopamine and are a clear "
            "risk signal, not a repurposing lead."
        ),
    ),
    Direction(
        "ACHE",
        protective=(),
        risk=(),
        confidence="ambiguous",
        rationale=(
            "Genuinely two-sided. Chlorpyrifos causes its toxicity by inhibiting AChE and "
            "is a PD risk factor, yet donepezil, rivastigmine and galantamine are "
            "inhibitors used for PD dementia, and rivastigmine reduces falls. The dose, "
            "duration and cell type differ from the toxin exposure. Do not let a score "
            "resolve this -- read the trial evidence per drug."
        ),
    ),
    Direction(
        "BCHE",
        protective=(),
        risk=(),
        confidence="ambiguous",
        rationale="Co-annotated with ACHE for the cholinesterase drugs; same two-sided argument.",
    ),
    Direction(
        "ESR1",
        protective=POSITIVE,
        risk=NEGATIVE,
        confidence="ambiguous",
        rationale=(
            "Oestrogen signalling is neuroprotective in toxin models and PD incidence is "
            "lower in women, but hormone-therapy trials are confounded and the systemic "
            "risk profile is substantial. Direction is arguable; certainty is not."
        ),
    ),
    Direction(
        "FTH1",
        protective=(),
        risk=(),
        confidence="ambiguous",
        rationale=(
            "Nigral iron accumulation is well established and iron chelation (deferiprone) "
            "has been trialled, but ferritin sequesters iron protectively -- so whether to "
            "raise or lower FTH1 itself does not follow from the chelation rationale."
        ),
    ),
    Direction(
        "KDR",
        protective=(),
        risk=(),
        confidence="ambiguous",
        rationale=(
            "VEGF signalling is neurotrophic in some PD models but drives blood-brain "
            "barrier leakage in others. The approved drugs here are oncology kinase "
            "inhibitors whose relevance is doubtful either way."
        ),
    ),
    Direction(
        "HPGDS",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="plausible",
        rationale=(
            "Haematopoietic prostaglandin D synthase drives PGD2-mediated neuroinflammation "
            "and reactive gliosis; inhibition is protective in models. Weak evidence base."
        ),
    ),
    Direction(
        "CD4",
        protective=NEGATIVE,
        risk=POSITIVE,
        confidence="ambiguous",
        rationale=(
            "T-cell infiltration contributes to nigral degeneration, but broad "
            "immunosuppression carries its own risk and regulatory T cells appear "
            "protective. Suppression is not straightforwardly beneficial."
        ),
    ),
    Direction(
        "CD2",
        protective=(),
        risk=(),
        confidence="unknown",
        rationale="Immune adhesion molecule; no PD-specific directional rationale.",
    ),
    # ADME and non-therapeutic annotations. Engagement here says nothing about
    # disease modification, so they are scored neutral rather than as support.
    *(
        Direction(
            gene,
            protective=(),
            risk=(),
            confidence="unknown",
            rationale=(
                "Metabolic, transport or carrier annotation rather than a therapeutic "
                "target; engaging it carries no directional PD rationale."
            ),
        )
        for gene in ("ALB", "ALAD", "CYP1B1", "CYP2E1", "GGT1", "TCEA1")
    ),
)

COLUMNS = ("gene", "protective_actions", "risk_actions", "confidence", "rationale")


def rows() -> list[tuple[str, str, str, str, str]]:
    """The table as loadable rows; action lists are joined with ``; ``."""
    return [
        (d.gene, "; ".join(d.protective), "; ".join(d.risk), d.confidence, d.rationale)
        for d in DIRECTIONS
    ]


def by_gene() -> dict[str, Direction]:
    return {d.gene: d for d in DIRECTIONS}


def classify(gene: str, action: str | None) -> str:
    """Label one (target, action_type) pair.

    Returns ``protective``, ``risk``, ``ambiguous`` or ``unknown``. Anything not
    explicitly curated is ``unknown`` -- absence of an entry is never support.
    """
    direction = by_gene().get(gene)
    if direction is None or not action:
        return "unknown"
    if action in direction.protective:
        return "protective"
    if action in direction.risk:
        return "risk"
    return "ambiguous" if direction.confidence == "ambiguous" else "unknown"
