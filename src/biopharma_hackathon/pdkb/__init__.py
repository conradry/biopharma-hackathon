"""A queryable evidence base for proposing repurposed drugs for Parkinson's.

Joins six layers into one DuckDB file: the PrimeKG disease neighbourhood, a
pathway enrichment, an environmental-toxin overlay, ChEMBL approved drugs with
their mechanism of action, those drugs' clinical trials, and the DrugCLIP
virtual screen. See :func:`build_database`.

Two things it is built to keep straight. Direction: a toxin and a drug can hit
one target in opposite ways, so :mod:`.mechanism` records which way opposes the
insult and the ranking gates on it. Provenance: every candidate can be unfolded
back into the evidence that produced it with :func:`evidence_graph`.
"""

from biopharma_hackathon.pdkb.build import DEFAULT_DB, build_database, connect
from biopharma_hackathon.pdkb.graph import EvidenceGraph, evidence_graph
from biopharma_hackathon.pdkb.mechanism import DIRECTIONS, Direction, classify

__all__ = [
    "DEFAULT_DB",
    "DIRECTIONS",
    "Direction",
    "EvidenceGraph",
    "build_database",
    "classify",
    "connect",
    "evidence_graph",
]
