"""Tests for the buried-surface layer.

Two things are checked that a contact-distance table cannot check itself:
that the burial measurements agree with the distances already committed, and
that burial and proximity are treated as different quantities.
"""

from __future__ import annotations

import pytest

from fcatlas.interface import (
    BURIED_AREA_THRESHOLD,
    burial_by_position,
    engineering_against_burial,
    load_interface,
    method_note,
)
from fcatlas.structure import load_complexes
from fcatlas.variants import load_variants

INTERFACES = load_interface()
COMPLEXES = load_complexes()
VARIANTS = {v.id: v for v in load_variants()}


def test_every_complex_measured():
    assert set(INTERFACES) == set(COMPLEXES)


@pytest.mark.parametrize("pdb", sorted(INTERFACES))
def test_contacts_agree_with_committed_distances(pdb):
    """The same positions, the same distances, from an independent pass."""
    detail = INTERFACES[pdb]
    committed = COMPLEXES[pdb].contacts
    assert set(detail.contact_positions) == set(committed)
    for eu, distance in committed.items():
        # the committed table takes the closest approach over both heavy chains
        assert detail.positions[eu].min_distance_any_chain == pytest.approx(
            distance, abs=0.01
        )


@pytest.mark.parametrize("pdb", sorted(INTERFACES))
def test_chains_and_genotype_match_the_deposition(pdb):
    detail, complex_ = INTERFACES[pdb], COMPLEXES[pdb]
    assert detail.fc_chains == complex_.fc_chains
    assert detail.partner_chains == complex_.partner_chains
    assert detail.fc_glycan_chains == complex_.fc_glycan_chains
    assert detail.deposited_fc_mutation == complex_.deposited_fc_mutation
    assert detail.deposited_partner_mutation == complex_.deposited_partner_mutation
    assert detail.resolution_angstrom == complex_.resolution_angstrom


@pytest.mark.parametrize("pdb", sorted(INTERFACES))
def test_areas_are_physically_sensible(pdb):
    detail = INTERFACES[pdb]
    for eu, position in detail.positions.items():
        assert position.sasa_free >= position.sasa_bound - 0.05, eu
        assert position.delta_sasa >= -0.05, eu
        assert -0.01 <= position.buried_fraction <= 1.01, eu
        assert position.delta_sasa == pytest.approx(
            position.sasa_free - position.sasa_bound, abs=0.15
        ), eu
    # A protein-protein interface buries a few hundred square angstroms a side.
    assert 400 < detail.total_buried_area < 1400


@pytest.mark.parametrize("pdb", sorted(INTERFACES))
def test_per_chain_maximum_is_the_reported_value(pdb):
    """The reported chain is the one that buries most, not an average."""
    for position in INTERFACES[pdb].positions.values():
        if not position.per_chain:
            continue
        best = max(position.per_chain.values(), key=lambda c: c["delta_sasa"])
        assert position.delta_sasa == pytest.approx(best["delta_sasa"], abs=0.05)


def test_representative_chain_distance_can_differ_from_the_closest_chain():
    """Distance and area are reported for the same copy of the residue.

    At EU 235 in the FcgammaRIIIb complex one heavy chain comes closer while
    the other buries more surface. Reporting the closest distance next to the
    largest area would describe two different residues in one row.
    """
    position = INTERFACES["1E4K"].positions[235]
    assert position.min_distance_any_chain <= position.min_distance
    assert position.min_distance_any_chain < position.min_distance
    assert position.delta_sasa == pytest.approx(
        position.per_chain[position.chain]["delta_sasa"], abs=0.05
    )


def test_burial_and_proximity_are_different_quantities():
    """EU 234 is closer to the receptor than EU 329 and buries far less.

    If this ever fails, either the measurement changed or the atlas is
    reporting distance dressed up as burial.
    """
    fcgr = INTERFACES["1E4K"]
    near, deep = fcgr.positions[234], fcgr.positions[329]
    assert near.min_distance < deep.min_distance
    assert deep.delta_sasa > 2 * near.delta_sasa
    assert deep.buried_fraction > near.buried_fraction


def test_asymmetric_positions_are_flagged():
    """Receptor binding to the Fc dimer is one-sided at some positions."""
    fcgr = INTERFACES["1E4K"]
    assert fcgr.positions[234].is_asymmetric
    assert fcgr.positions[329].is_asymmetric
    # and the asymmetry points in opposite directions
    assert fcgr.positions[234].chain != fcgr.positions[329].chain


def test_yte_third_substitution_makes_no_contact():
    """T256E lies outside the contact shell of the structure that carries it."""
    fcrn = INTERFACES["4N0U"]
    assert fcrn.deposited_fc_mutation is not None
    assert "T256E" in fcrn.deposited_fc_mutation.replace(" ", "")
    position = fcrn.positions[256]
    assert position.min_distance is None
    assert not position.in_contact
    assert position.delta_sasa < 20
    # the other two do touch it
    assert fcrn.positions[252].in_contact
    assert fcrn.positions[254].in_contact


def test_glycan_is_part_of_the_interface():
    for pdb in ("1E4K", "5XJE"):
        assert INTERFACES[pdb].glycan_buried_area > 20


def test_engineered_partner_is_recorded():
    assert INTERFACES["5XJE"].deposited_partner_mutation is not None
    assert not INTERFACES["5XJE"].wild_type_partner
    assert INTERFACES["1E4K"].wild_type_partner


def test_variant_lookup_finds_its_own_positions():
    lala = VARIANTS["LALA"]
    rows = INTERFACES["1E4K"].for_variant(lala)
    assert {p.eu for p in rows} == {234, 235}
    assert rows[0].delta_sasa >= rows[-1].delta_sasa
    assert INTERFACES["1E4K"].variant_buried_area(lala) > 100


def test_burial_by_position_reports_the_largest_and_says_where():
    table = burial_by_position(INTERFACES)
    for eu, row in table.items():
        assert row["pdb"] in INTERFACES
        assert row["delta_sasa"] == pytest.approx(
            INTERFACES[row["pdb"]].positions[eu].delta_sasa
        )


def test_most_engineered_position_is_not_most_buried():
    """The field's attention and the interface's geometry are not aligned."""
    rows = engineering_against_burial(load_variants(), INTERFACES)
    measured = [r for r in rows if r["delta_sasa"]]
    most_engineered = max(rows, key=lambda r: r["records"])
    most_buried = max(measured, key=lambda r: r["delta_sasa"])
    assert most_engineered["eu"] != most_buried["eu"]
    # and the most buried FcRn position is one nobody in this set engineers
    fcrn_top = max(
        (r for r in measured if r["pdb"] == "4N0U"), key=lambda r: r["delta_sasa"]
    )
    assert fcrn_top["records"] == 0


def test_threshold_is_stated_not_implied():
    assert BURIED_AREA_THRESHOLD == 10.0
    note = method_note()
    assert note["probe_radius_angstrom"] == 1.4
    assert "shrake" in note["sasa_algorithm"].lower()
    assert note["contact_cutoff_angstrom"] == 5.0
