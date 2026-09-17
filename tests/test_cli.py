"""The command line is the part most people touch, so it is tested as text.

These tests assert on what a user reads, not on internal calls: the exit code,
the presence of the numbers that matter, and the presence of the caveats that
stop a number from being read as more than it is. A silent caveat is a bug here,
which is why the provenance and burial commands are checked for their scope lines
as well as their values.
"""

from __future__ import annotations

import re

import pytest

from fcatlas.cli import main


def run(capsys, *argv):
    code = main(list(argv))
    out = capsys.readouterr().out
    return code, out


# ------------------------------------------------------------------ burial
def test_burial_ranks_by_area(capsys):
    code, out = run(capsys, "burial", "--pdb", "1E4K", "--top", "5")
    assert code == 0
    assert "1E4K" in out
    # EU329 buries the most surface in 1E4K and must therefore lead the table.
    rows = re.findall(r"EU(\d+) \w+ chain \w+: [\d.]+ A, buries ([\d.]+) A", out)
    assert rows, out
    assert rows[0][0] == "329"
    assert rows[0][1] == "120.5"
    areas = [float(a) for _, a in rows]
    assert areas == sorted(areas, reverse=True)


def test_burial_reports_asymmetry(capsys):
    _, out = run(capsys, "burial", "--pdb", "1E4K", "--top", "3")
    # EU329 buries 120.5 A^2 on one heavy chain and 5.2 on the other. A reader
    # who does not see both numbers will average two that should not be averaged,
    # so the per-chain split has to be printed, not just alluded to in the note.
    assert "chain A 120.5 A^2" in out
    assert "chain B 5.2 A^2" in out


def test_burial_against_density_shows_both_extremes(capsys):
    code, out = run(capsys, "burial", "--against-density")
    assert code == 0
    assert "253" in out and "121.2" in out  # buried, never engineered
    assert "234" in out  # engineered, barely buried
    assert "records" in out.lower()


def test_burial_states_its_scope(capsys):
    _, out = run(capsys, "burial", "--pdb", "4N0U")
    low = out.lower()
    assert "shrake" in low or "1.40" in out
    assert "4N0U" in out


# ------------------------------------------------------------------ partners
def test_partners_names_the_contacting_residue(capsys):
    code, out = run(capsys, "partners", "LALA-PG", "--pdb", "1E4K")
    assert code == 0
    assert "TRP" in out.upper()
    assert "3.26" in out
    assert "329" in out


def test_partners_says_so_when_nothing_is_touched(capsys):
    # A CH3 heterodimerization variant touches no receptor. An empty table
    # would look identical to a failed lookup, so the command must use words.
    code, out = run(capsys, "partners", "knobs-into-holes", "--pdb", "1E4K")
    assert code == 0
    assert out.strip()
    low = out.lower()
    assert "no " in low or "none" in low


def test_partners_rejects_an_unknown_variant(capsys):
    code, _ = run(capsys, "partners", "not-a-real-variant", "--pdb", "1E4K")
    assert code != 0


# ------------------------------------------------------------------ provenance
def test_provenance_prints_both_genotypes(capsys):
    code, out = run(capsys, "provenance")
    assert code == 0
    # The receptor in 5XJE is engineered at four positions. That has to be
    # visible, because it is the only FcgammaRIIIa complex in the set.
    assert "5XJE" in out
    assert "N56Q" in out or "N187Q" in out
    # And the Fc in 4N0U is not wild type.
    assert "M252Y" in out


def test_provenance_identifies_1t89_as_iiib(capsys):
    _, out = run(capsys, "provenance")
    block = [ln for ln in out.splitlines() if "1T89" in ln]
    assert block
    joined = " ".join(block).lower()
    assert "iiib" in joined or "iiib" in out.lower()


def test_provenance_reports_resolution(capsys):
    _, out = run(capsys, "provenance")
    for value in ("3.2", "3.5", "2.4", "3.8"):
        assert value in out


# ------------------------------------------------------------------ interface
def test_interface_carries_area_next_to_distance(capsys):
    code, out = run(capsys, "interface", "LALA-PG")
    assert code == 0
    assert "3.26" in out or "120.5" in out


@pytest.mark.parametrize(
    "argv",
    [
        ("list",),
        ("hotspots",),
        ("structures",),
        ("provenance",),
        ("burial",),
        ("validate",),
    ],
)
def test_commands_exit_clean(capsys, argv):
    code, out = run(capsys, *argv)
    assert code == 0
    assert out.strip()


def test_partners_keeps_a_position_that_buries_but_touches_nothing(capsys):
    # T256E in 4N0U loses 15.1 A^2 and has no partner atom inside the cutoff.
    # Reporting only contacting positions would make YTE look like a two-position
    # variant at this interface.
    code, out = run(capsys, "partners", "YTE", "--pdb", "4N0U")
    assert code == 0
    assert "EU256" in out
    assert "15.1" in out
    assert "no partner atom" in out
