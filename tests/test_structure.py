"""Structural claims are distances, and their provenance is part of the claim."""

from fcatlas import (
    classify,
    fcrn_provenance,
    find,
    interface_report,
    load_complexes,
    load_variants,
    pymol_script,
)

CX = load_complexes()


def test_expected_complexes_present():
    assert sorted(CX) == ["1E4K", "1T89", "4N0U", "5XJE"]


def test_contact_counts_match_measured_values():
    assert len(CX["1E4K"].contacts) == 19
    assert len(CX["1T89"].contacts) == 18
    assert len(CX["5XJE"].contacts) == 23
    assert len(CX["4N0U"].contacts) == 16


def test_cutoff_is_uniform_and_stated():
    assert {c.cutoff for c in CX.values()} == {5.0}


def test_lower_hinge_dominates_the_fcgr_interface():
    for pdb in ("1E4K", "1T89", "5XJE"):
        assert CX[pdb].at_interface(235)
    # 1T89 begins at 235; 1E4K resolves 234 as well.
    assert CX["1E4K"].at_interface(234)
    assert not CX["1T89"].at_interface(234)


def test_closest_fcrn_contacts_are_310_and_434():
    c = CX["4N0U"]
    closest = sorted(c.contacts.items(), key=lambda kv: kv[1])[:2]
    assert {eu for eu, _ in closest} == {310, 434}
    assert all(d < 2.1 for _, d in closest)


def test_no_deposited_fc_fcrn_complex_carries_a_wild_type_fc():
    """The finding that limits every structural statement about FcRn here."""
    fcrn = [c for c in CX.values() if "FcRn" in c.partner]
    assert fcrn and all(not c.is_wild_type for c in fcrn)
    assert "YTE" in CX["4N0U"].fc_genotype
    assert CX["4N0U"].caveat


def test_fcrn_provenance_survey_is_carried_with_the_data():
    s = fcrn_provenance()
    entries = s["entries"]
    counts = s["_counts"]
    assert counts["entries_surveyed"] == len(entries) == 21
    assert counts["fcrn_containing"] == 9
    assert counts["fc_fcrn_complexes"] == 4
    assert counts["fc_fcrn_complexes_with_wild_type_fc"] == 0
    assert every_category_is_known(entries)


def every_category_is_known(entries):
    allowed = {"fc_fcrn_complex", "fcrn_without_fc", "unbound_fc"}
    return all(e["category"] in allowed for e in entries)


def test_the_four_fc_fcrn_complexes_are_all_engineered():
    entries = fcrn_provenance()["entries"]
    cx = {e["pdb"]: e for e in entries if e["category"] == "fc_fcrn_complex"}
    assert sorted(cx) == ["4N0U", "6WNA", "6WOL", "7Q15"]
    for pdb, e in cx.items():
        assert e["has_igg_fc"]
        assert e["fc_chains_engineered"], f"{pdb} should have engineered Fc chains"
        assert e["fc_chains_wild_type"] == [], f"{pdb} unexpectedly has wild-type Fc"


def test_wild_type_human_fc_appears_only_unbound():
    entries = fcrn_provenance()["entries"]
    wt = {e["pdb"]: e for e in entries if e["fc_chains_wild_type"]}
    assert set(wt) == {"4WI2", "6BZ4", "1L6X"}
    assert all(e["category"] == "unbound_fc" for e in wt.values())


def test_fcrn_entries_without_fc_are_albumin_or_viral():
    entries = fcrn_provenance()["entries"]
    no_fc = [e for e in entries if e["category"] == "fcrn_without_fc"]
    assert sorted(e["pdb"] for e in no_fc) == ["4K71", "4N0F", "6QIO", "6QIP", "9NAV"]
    for e in no_fc:
        assert not e["has_igg_fc"]
        assert "albumin" in e["composition"] or "capsid" in e["composition"]


def test_alanine_scan_series_is_unbound_fc():
    by_pdb = {e["pdb"]: e for e in fcrn_provenance()["entries"]}
    assert by_pdb["4WI3"]["fc"][0][1] == ["I253A"]
    for i in range(2, 10):
        assert by_pdb[f"4WI{i}"]["category"] == "unbound_fc"
    assert "inferred from engineered Fc" in fcrn_provenance()["_finding"]


def test_glycosylated_model_resolves_more_interface_than_deglycosylated():
    assert len(CX["5XJE"].contacts) > len(CX["1T89"].contacts)
    added = set(CX["5XJE"].contacts) - set(CX["1T89"].contacts)
    assert {268, 294, 295, 296}.issubset(added)


def test_interface_report_reports_absence_explicitly():
    v = find("E430G")[0]
    r = interface_report(v)
    assert all(
        all(d is None for d in row["distances"].values()) for row in r.rows
    )
    assert "no measured interface contact" in r.summary()


def test_classification_counts():
    rows = classify()
    assert len(rows) == len(load_variants())
    fcgr = [r for r in rows if r["fcgr_interface"]]
    fcrn = [r for r in rows if r["fcrn_interface"]]
    neither = [r for r in rows if not r["fcgr_interface"] and not r["fcrn_interface"]]
    assert len(fcgr) == 42
    assert len(fcrn) == 12
    assert len(neither) == 24


def test_heterodimer_variants_sit_outside_both_interfaces():
    rows = {r["id"]: r for r in classify()}
    for vid in ("KiH-orig-knob", "DuoBody", "DuoBody-partner", "DEKK-DE", "DEKK-KK"):
        assert vid in rows
        assert not rows[vid]["fcgr_interface"]
        assert not rows[vid]["fcrn_interface"]


def test_pymol_script_is_selfcontained_and_names_unresolved_residues():
    v = find("LALA-PG")[0]
    s = pymol_script(v, pdb="4N0U")
    assert "fetch 4N0U" in s
    assert "LALA-PG" in s
    # None of 234/235/329 is resolved at the FcRn interface, and the script
    # must say so rather than silently drawing nothing.
    assert "none of this variant's positions contact the partner" in s
    assert "CAVEAT" in s and "not wild type" in s


def test_pymol_script_is_deterministic():
    v = find("YTE")[0]
    assert pymol_script(v, pdb="4N0U") == pymol_script(v, pdb="4N0U")
