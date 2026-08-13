"""Tests for the toxin-target parser and its hand-curated EC mapping.

EC_MAP is curation, not code -- an edit there changes results silently, and the
mechanistic anchors (MPTP reaching monoamine oxidase, rotenone reaching complex
I) are exactly what a careless edit would drop. The fixture reproduces the
awkward cases in the real file: a protein listed twice under both a symbol and
an EC, a family EC that expands to several genes, non-human and ambiguous EC
classes, and a mangled symbol.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_toxin_targets import EC_MAP, parse_toxins

HEADER = (
    "Toxin,PubMedMentions,CID,HasHumanData,MatchedSections,"
    "NumGeneTargets,TopGeneTargets,SampleExcerpt\n"
)


def _write(tmp_path: Path, rows: list[tuple[str, str]]) -> Path:
    path = tmp_path / "toxins.csv"
    body = "".join(f'{name},5,123,True,,25,"{targets}",\n' for name, targets in rows)
    path.write_text(HEADER + body)
    return path


def _entries(path: Path, toxin: str) -> dict[str, float]:
    """Effective evidence per gene, the way the build aggregates it."""
    (parsed,) = [t for t in parse_toxins(path) if t.name == toxin]
    totals: dict[str, float] = {}
    for entry in parsed.entries:
        totals[entry.symbol] = totals.get(entry.symbol, 0.0) + entry.count / entry.expansion
    return totals


def test_mptp_reaches_monoamine_oxidase(tmp_path: Path) -> None:
    """The MPTP anchor: MAO is named only by EC, and it is the causal target."""
    path = _write(tmp_path, [("mptp", "monoamine oxidase [EC:1.4.3.4] (481)")])
    assert _entries(path, "mptp") == {"MAOA": 240.5, "MAOB": 240.5}


def test_rotenone_reaches_complex_i(tmp_path: Path) -> None:
    path = _write(tmp_path, [("rotenone", "nadh dehydrogenase [EC:1.6.99.3] (463)")])
    genes = _entries(path, "rotenone")
    assert set(genes) == {f"NDUFS{n}" for n in (1, 2, 3, 7, 8)} | {"NDUFV1", "NDUFV2"}
    assert sum(genes.values()) == pytest.approx(463.0)


def test_symbol_and_ec_for_one_protein_are_summed_not_duplicated(tmp_path: Path) -> None:
    """Lead lists catalase twice; the counts add, and the gene appears once."""
    path = _write(tmp_path, [("lead", "catalase [CAT] (576); catalase [EC:1.11.1.6] (234)")])
    assert _entries(path, "lead") == {"CAT": 810.0}


def test_family_evidence_is_split_not_replicated(tmp_path: Path) -> None:
    """A family EC must not outvote a specific symbol by sheer expansion."""
    path = _write(
        tmp_path,
        [
            (
                "paraquat",
                "superoxide dismutase [EC:1.15.1.1] (969); superoxide dismutase 2 [SOD2] (188)",
            )
        ],
    )
    genes = _entries(path, "paraquat")
    assert genes == {"SOD1": 323.0, "SOD2": 511.0, "SOD3": 323.0}
    assert sum(genes.values()) == pytest.approx(969 + 188)


def test_non_human_and_ambiguous_classes_are_dropped_with_a_reason(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        [("manganese", "photosystem ii [EC:1.10.3.9] (1148); peroxidase [EC:1.11.1.7] (845)")],
    )
    (parsed,) = parse_toxins(path)
    assert parsed.entries == []
    assert sorted(reason for _, reason in parsed.unmapped) == ["ambiguous", "non_human"]


def test_mangled_symbol_is_corrected(tmp_path: Path) -> None:
    path = _write(tmp_path, [("dieldrin", "aefer(h) [FER1HCH] (458)")])
    assert _entries(path, "dieldrin") == {"FTH1": 458.0}


def test_listed_count_covers_dropped_entries(tmp_path: Path) -> None:
    """The mapped-fraction denominator must include what was thrown away."""
    path = _write(
        tmp_path,
        [
            (
                "mercury",
                "catalase [CAT] (363); mercuric reductase [MERA] (211); "
                "photosystem ii [EC:1.10.3.9] (100)",
            )
        ],
    )
    (parsed,) = parse_toxins(path)
    assert parsed.listed_count == 363 + 211 + 100


def test_rows_without_targets_are_skipped(tmp_path: Path) -> None:
    path = tmp_path / "toxins.csv"
    path.write_text(HEADER + "air pollution,12,,False,,0,,\n")
    assert parse_toxins(path) == []


def test_every_ec_class_declares_a_confidence() -> None:
    valid = {"exact", "family", "complex", "ambiguous", "non_human"}
    assert all(confidence in valid for _, confidence in EC_MAP.values())
    for symbols, confidence in EC_MAP.values():
        assert bool(symbols) == (confidence in {"exact", "family", "complex"})
