"""The dataset must never disagree with the reference sequence in silence."""

import pytest

from fcatlas import (
    INTENT_LABELS,
    by_intent,
    find,
    load_complexes,
    load_variants,
    numbering,
    position_index,
    validate_all,
)

VARIANTS = load_variants()


def test_dataset_size():
    assert len(VARIANTS) == 74


def test_no_unexplained_validation_problems():
    """Every wild-type residue implied by a mutation string is checked.

    An empty result is the contract: a substitution written L234A must find a
    leucine at EU 234 in the parent isotype, or the record is wrong.
    """
    assert validate_all() == []


def test_every_record_has_a_source():
    missing = [v.id for v in VARIANTS if not v.source.strip()]
    assert missing == []


def test_ids_are_unique():
    ids = [v.id for v in VARIANTS]
    assert len(set(ids)) == len(ids)


def test_intents_are_declared():
    unknown = sorted({v.intent for v in VARIANTS} - set(INTENT_LABELS))
    assert unknown == []


def test_substitution_count():
    n = sum(len(v.substitutions) for v in VARIANTS)
    assert n == 174


def test_disputed_records_are_flagged_and_explained():
    disputed = [v for v in VARIANTS if v.is_disputed]
    assert len(disputed) == 2
    for v in disputed:
        assert v.notes, f"{v.id} carries a disputed substitution with no explanation"


def test_igg2_sigma_dispute_is_the_numbering_collision():
    v = find("IgG2-sigma")[0]
    assert any(s.eu == 234 for s in v.disputed)
    assert numbering("IgG2").residue(234) is None


def test_eculizumab_hybrid_is_curated_as_igg4():
    """Its implied wild-type residues are the IgG4 hinge residues.

    IMGT labels this class of hybrid G2G4v1, which invites reading it in the
    IgG2 frame; validating against sequence forces the IgG4 frame instead.
    """
    v = find("eculizumab-hybrid")[0]
    assert v.isotype == "IgG4"
    assert v.validate() == []


def test_charge_pair_kk_dispute_is_allotypic():
    v = find("charge-pair-KK")[0]
    assert any(s.eu == 356 for s in v.disputed)
    assert "allotype" in v.notes.lower()


def test_position_index_finds_the_dominant_hotspot():
    idx = position_index()
    ranked = sorted(idx.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    top_eu, top_variants = ranked[0]
    assert top_eu == 235
    assert len(top_variants) == 20
    assert len(ranked[1][1]) == 17 and ranked[1][0] == 234


def test_distinct_positions_touched():
    assert len(position_index()) == 60


def test_by_intent_partitions_the_dataset():
    groups = by_intent()
    assert sum(len(g) for g in groups.values()) == len(VARIANTS)
    assert len(groups["silence"]) == 27


def test_mutated_sequence_applies_only_validated_substitutions():
    v = find("LALA")[0]
    wt = numbering("IgG1").sequence
    mut = v.mutated_sequence()
    assert len(mut) == len(wt)
    diffs = [i for i, (a, b) in enumerate(zip(wt, mut)) if a != b]
    assert len(diffs) == 2
    n = numbering("IgG1")
    assert {n.positions()[i] for i in diffs} == {234, 235}


@pytest.mark.parametrize("query", ["LALA-PG", "lala-pg", "P329G", "YTE"])
def test_find_is_case_insensitive_and_matches_mutations(query):
    assert find(query)


def test_every_source_is_a_resolvable_style_reference():
    """Sources must be DOI URLs, or an explicitly allowed non-DOI reference.

    Every DOI in this file was checked against Crossref during curation to
    confirm that it resolves to the paper the record claims. This test cannot
    re-check resolution offline, but it does keep the format enforceable and
    keeps the allow-list of non-DOI sources visible and small.
    """
    allowed_non_doi = {
        "https://patents.google.com/patent/WO2015195498A1/en",
    }
    for variant in load_variants():
        source = variant.source
        assert source, f"{variant.id} has no source"
        if source in allowed_non_doi:
            continue
        assert source.startswith("https://doi.org/10."), (
            f"{variant.id} source is neither a DOI URL nor allow-listed: {source}"
        )


def test_structure_dois_are_doi_urls():
    for pdb_id, complex_ in load_complexes().items():
        assert complex_.doi.startswith("https://doi.org/10."), pdb_id
