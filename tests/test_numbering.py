"""EU numbering is the load-bearing assumption of this package.

If these tests fail, every downstream position claim is suspect.
"""

import pytest

from fcatlas import EU_START, ISOTYPES, align_isotypes, constant_region, domain_of, numbering

# Landmarks with a residue identity fixed by the literature or by the
# UniProt reference sequence itself. Positions are EU; residues are IgG1.
IGG1_LANDMARKS = {
    118: "A",   # first residue of the secreted constant region
    220: "C",   # light-chain disulfide cysteine
    226: "C",   # hinge cysteine
    229: "C",   # hinge cysteine
    234: "L",   # LALA
    235: "L",   # LALA
    236: "G",   # GRLR / LALAGA
    237: "G",
    239: "S",   # S239D
    252: "M",   # YTE
    254: "S",
    256: "T",
    265: "D",   # D265A
    267: "S",
    270: "D",   # D270A
    297: "N",   # N297 glycosylation sequon
    298: "S",
    299: "T",
    322: "K",   # K322A
    326: "K",   # complement enhancement
    328: "L",   # LALA-LR
    329: "P",   # P329G proline sandwich
    330: "A",
    331: "P",   # P331S
    332: "I",   # I332E
    345: "E",   # E345K hexamerization
    366: "T",   # knob-into-hole T366
    368: "L",
    407: "Y",   # hole Y407V
    409: "K",   # IgG4 R409K counterpart
    428: "M",   # LS / MST-HN
    429: "H",
    434: "N",   # N434S / N434A
    435: "H",
    447: "K",   # C-terminal lysine
}


def test_eu_start_is_118():
    assert EU_START == 118


def test_igg1_landmark_residues():
    n = numbering("IgG1")
    wrong = {eu: n.residue(eu) for eu, aa in IGG1_LANDMARKS.items() if n.residue(eu) != aa}
    assert wrong == {}


def test_igg1_range_and_length():
    n = numbering("IgG1")
    pos = n.positions()
    assert pos[0] == 118
    assert pos[-1] == 447
    assert len(constant_region("IgG1")) == 330


def test_all_isotypes_numbered():
    lengths = {i: len(constant_region(i)) for i in ISOTYPES}
    assert lengths == {"IgG1": 330, "IgG2": 326, "IgG3": 377, "IgG4": 327}


def test_igg3_has_unnumbered_hinge_insertion():
    # IgG3 carries the extended hinge; those residues have no EU counterpart
    # in the IgG1 frame and must be reported as unnumbered, not silently
    # renumbered.
    assert len(numbering("IgG3").unnumbered_indices) == 47
    assert numbering("IgG1").unnumbered_indices == []


def test_igg2_has_no_residue_at_eu234():
    """The finding this package exists to make unavoidable.

    Aligned to IgG1, IgG2 has a gap at EU 234; its valine sits at EU 235.
    Primary literature writes the IgG2-sigma substitution as V234A while
    IMGT writes V235A. Both describe the same physical residue.
    """
    igg2 = numbering("IgG2")
    assert igg2.residue(234) is None
    assert igg2.residue(235) == "V"
    # IgG1 and IgG4, by contrast, both have residues at 234.
    assert numbering("IgG1").residue(234) == "L"
    assert numbering("IgG4").residue(234) == "F"


def test_cross_isotype_alignment_preserves_landmarks():
    aln = align_isotypes()
    assert len(aln) == 330
    for eu, aa in IGG1_LANDMARKS.items():
        cols = [col for pos, col in aln if pos == eu]
        assert len(cols) == 1, f"EU {eu} not represented exactly once"
        assert cols[0]["IgG1"] == aa


def test_eu358_is_allotype_dependent_in_reference():
    """EU 356/358 are the classical IgG1 allotype pair.

    The reference sequence used here is G1m(f) (D356/L358). Any variant
    written against G1m(1) (E356/M358) will disagree, which is why those
    records are carried as disputed rather than corrected.
    """
    n = numbering("IgG1")
    assert n.residue(356) == "D"
    assert n.residue(358) == "L"


@pytest.mark.parametrize(
    "eu,expected",
    [(118, "CH1"), (230, "hinge"), (300, "CH2"), (400, "CH3")],
)
def test_domain_assignment(eu, expected):
    assert domain_of(eu) == expected


def test_numbering_is_cached_and_immutable_across_calls():
    a, b = numbering("IgG1"), numbering("IgG1")
    assert a.residue(297) == b.residue(297) == "N"
