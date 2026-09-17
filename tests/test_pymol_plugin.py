"""Tests for the PyMOL plugin, without PyMOL.

A structural plugin is easy to ship broken: the selection strings are strings,
so a wrong chain letter or a stale residue list draws a confident picture of
the wrong atoms and nothing raises. These tests install a fake ``pymol``
module that records every call, then assert on what the plugin asked for.

The assertions are about correctness of the claim, not cosmetics: that the
chains coloured as receptor are the receptor chains in that deposition, that
the residues marked as a variant are the residues that variant changes, that
burial written into the B-factor column matches the committed measurement,
and that sugars are drawn as spheres rather than lines, since their bond
connectivity is inferred by the viewer rather than deposited.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pymol_plugin"))


class FakeCmd:
    """Records calls; answers the few queries the plugin makes."""

    def __init__(self):
        self.calls: list[tuple] = []
        self.objects: list[str] = []
        self.selections: dict[str, str] = {}
        self.altered: list[tuple[str, str]] = []
        self.atom_counts: dict[str, int] = {}
        self.extended: dict[str, object] = {}
        self.auto_arg = [{}, {}, {}]

    def _record(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))

    def __getattr__(self, name):
        def call(*args, **kwargs):
            self._record(name, *args, **kwargs)

        return call

    # queries and side effects the plugin depends on
    def fetch(self, code, name=None, **kwargs):
        self._record("fetch", code, name, **kwargs)
        self.objects.append(name or code)

    def get_names(self, kind="objects", *args, **kwargs):
        self._record("get_names", kind)
        return list(self.objects)

    def select(self, name, selection, **kwargs):
        self._record("select", name, selection)
        self.selections[name] = selection
        return self.atom_counts.get(name, 1)

    def count_atoms(self, selection, **kwargs):
        self._record("count_atoms", selection)
        return self.atom_counts.get(selection, 1)

    def alter(self, selection, expression, **kwargs):
        self._record("alter", selection, expression)
        self.altered.append((selection, expression))

    def delete(self, name):
        self._record("delete", name)
        if name == "all":
            self.objects.clear()
            self.selections.clear()

    def extend(self, name, function):
        self._record("extend", name)
        self.extended[name] = function

    def of(self, name):
        return [c for c in self.calls if c[0] == name]

    def kwargs_of(self, name, index=0):
        return self.of(name)[index][2]

    def args_of(self, name):
        return [c[1] for c in self.calls if c[0] == name]


@pytest.fixture()
def plugin(monkeypatch):
    """The plugin module, wired to a fresh recorder."""
    fake_cmd = FakeCmd()
    module = types.ModuleType("pymol")
    module.cmd = fake_cmd
    monkeypatch.setitem(sys.modules, "pymol", module)
    for stale in [m for m in sys.modules if m.startswith("fcatlas_pymol")]:
        monkeypatch.delitem(sys.modules, stale, raising=False)
    import fcatlas_pymol

    fcatlas_pymol._CACHE.clear()
    return fcatlas_pymol, fake_cmd


# --------------------------------------------------------------------- wiring


def test_data_is_found_and_complete(plugin):
    module, _ = plugin
    assert module.data_dir().is_dir()
    assert len(module.variants()) > 50
    assert set(module.structures()) >= {"1E4K", "1T89", "5XJE", "4N0U"}
    assert module.interface_detail()["complexes"]


def test_commands_register(plugin):
    module, cmd = plugin
    names = module.register_commands()
    assert "fcatlas_show" in names
    assert set(cmd.extended) == set(module.COMMANDS)
    for function in cmd.extended.values():
        assert callable(function)


def test_plugin_entry_point_runs_without_qt(plugin):
    module, cmd = plugin
    module.__init_plugin__(None)  # no pymol.plugins in the fake module
    assert "fcatlas_show" in cmd.extended


def test_variant_lookup_is_forgiving(plugin):
    module, _ = plugin
    assert module.find_variant("lala")["id"] == "LALA"
    assert module.find_variant("LALA-PG")["id"] == "LALA-PG"
    with pytest.raises(KeyError):
        module.find_variant("not-a-variant")


def test_half_life_variants_default_to_the_fcrn_complex(plugin):
    module, _ = plugin
    assert module.default_structure(module.find_variant("YTE")) == "4N0U"
    assert module.default_structure(module.find_variant("LALA")) == "1E4K"


# ---------------------------------------------------------------- scene claims


def test_show_colours_the_deposited_chains_not_guessed_ones(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    structure = module.structures()["1E4K"]
    # cmd.color(colour, selection)
    coloured = {colour: selection for colour, selection in cmd.args_of("color")}
    partner_selection = coloured[module.PARTNER_COLOR]
    for chain in structure["partner_chains"]:
        assert f"chain {chain}" in partner_selection or chain in partner_selection
    # the Fc chains must not be coloured as partner
    for chain in structure["fc_chains"]:
        assert f"chain {chain} " not in partner_selection.replace("+", " ") + " "


def test_1t89_receptor_chain_is_c_not_b(plugin):
    """The regression this catches: 1T89 has Fc on A and B, receptor on C."""
    module, cmd = plugin
    module.fcatlas_show("LALA", "1T89")
    assert module.structures()["1T89"]["partner_chains"] == ["C"]
    assert module.structures()["1T89"]["fc_chains"] == ["A", "B"]
    variant_selection = cmd.selections["fc_LALA_1T89_variant"]
    assert "chain A+B" in variant_selection


def test_variant_selection_holds_exactly_the_changed_positions(plugin):
    module, cmd = plugin
    record = module.find_variant("LALA-PG")
    module.fcatlas_show("LALA-PG", "1E4K")
    selection = cmd.selections["fc_LALA_PG_1E4K_variant"]
    listed = selection.split("resi ")[1].split()[0]
    assert sorted(int(x) for x in listed.split("+")) == sorted(
        module.positions_of(record)
    )


def test_interface_selection_matches_the_committed_contacts(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    selection = cmd.selections["fc_LALA_1E4K_interface"]
    listed = {int(x) for x in selection.split("resi ")[1].split()[0].split("+")}
    assert listed == {int(k) for k in module.structures()["1E4K"]["contacts"]}


def test_sugars_are_spheres_and_never_lines(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    shown = {args[0]: args[1] for args in cmd.args_of("show")}
    glycan_shows = [rep for rep, sel in shown.items() if "glycan" in str(sel)]
    assert "spheres" in glycan_shows
    assert not any(rep == "lines" for rep in shown)
    assert "D+E" in cmd.selections["fc_LALA_1E4K_glycans"]


def test_solvent_is_removed_before_measurement_is_shown(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    assert any("solvent" in args[0] for args in cmd.args_of("remove"))


def test_scene_fetches_the_requested_entry_synchronously(plugin):
    module, cmd = plugin
    module.fcatlas_show("YTE")
    code, name = cmd.of("fetch")[0][1][:2]
    assert code == "4N0U"
    assert name == "fc_YTE_4N0U"
    assert cmd.of("fetch")[0][2]["async_"] == 0


def test_scene_report_states_provenance(plugin):
    module, _ = plugin
    report = module.fcatlas_show("YTE", "4N0U")
    assert "M252Y" in report and "T256E" in report
    assert "3.8 A" in report
    assert "albumin" in report.lower()
    assert "P01857" in report or "genotype" in report


def test_engineered_receptor_is_declared(plugin):
    module, _ = plugin
    report = module.fcatlas_provenance("5XJE")
    assert "F176V" in report
    assert "FcgammaRIIIa" in report
    # and the wild-type receptor entries say so
    assert "wild type" in module.fcatlas_provenance("1E4K")


# -------------------------------------------------------------------- burial


def test_burial_writes_measured_areas_into_b_factors(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    module.fcatlas_burial("1E4K")
    zeroed = [a for a in cmd.altered if a[1] == "b=0.0"]
    assert zeroed, "the column must be cleared before it is used"
    written = {}
    for selection, expression in cmd.altered:
        if expression == "b=0.0":
            continue
        eu = int(selection.split("resi ")[1])
        written[eu] = float(expression.split("=")[1])
    detail = module.interface_detail()["complexes"]["1E4K"]["positions"]
    for eu, value in written.items():
        assert value == pytest.approx(detail[str(eu)]["delta_sasa"])
    assert written[329] > written[234]


def test_burial_writes_to_the_chain_that_buries_the_surface(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    module.fcatlas_burial("1E4K")
    detail = module.interface_detail()["complexes"]["1E4K"]["positions"]
    for selection, expression in cmd.altered:
        if expression == "b=0.0":
            continue
        chain = selection.split("chain ")[1].split()[0]
        eu = selection.split("resi ")[1]
        assert chain == detail[eu]["chain"]


def test_burial_colours_by_b_and_nothing_else(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    module.fcatlas_burial("1E4K")
    expression, palette = cmd.of("spectrum")[0][1][:2]
    assert expression == "b"
    assert palette == module.BURIAL_PALETTE
    assert cmd.of("spectrum")[0][2]["minimum"] == 0


def test_burial_without_a_loaded_object_says_so(plugin):
    module, _ = plugin
    assert "Nothing loaded" in module.fcatlas_burial("1E4K")


# ------------------------------------------------------------------- partners


def test_partners_selects_named_receptor_residues(plugin):
    module, cmd = plugin
    module.fcatlas_show("LALA", "1E4K")
    report = module.fcatlas_partners("LALA", "1E4K")
    selection = cmd.selections["fc_LALA_1E4K_partner_side"]
    assert "chain C" in selection
    detail = module.interface_detail()["complexes"]["1E4K"]["positions"]
    expected = {
        p["seq"] for eu in ("234", "235") for p in detail[eu]["partners"]
    }
    listed = {
        int(x)
        for group in selection.split("resi ")[1:]
        for x in group.split(")")[0].split("+")
    }
    assert listed == expected
    assert "HIS116" in report.replace(" ", "") or "HIS116" in report


def test_partners_reports_absence_rather_than_inventing_contacts(plugin):
    module, _ = plugin
    report = module.fcatlas_partners("M252E", "1E4K") if any(
        v["id"] == "M252E" for v in _plugin_variants(module)
    ) else module.fcatlas_partners("YTE", "1E4K")
    assert "no position in contact" in report or "no contact" in report


def _plugin_variants(module):
    return module.variants()


# ---------------------------------------------------------------- text reports


def test_interface_report_separates_distance_from_burial(plugin):
    module, _ = plugin
    report = module.fcatlas_interface("LALA-PG", "1E4K")
    assert "buried" in report.lower()
    assert "120.5" in report  # EU329
    assert "asymmetric" in report


def test_interface_report_names_positions_that_do_not_touch(plugin):
    module, _ = plugin
    report = module.fcatlas_interface("YTE", "4N0U")
    assert "no contact within the cutoff" in report


def test_hotspots_puts_records_next_to_burial(plugin):
    module, _ = plugin
    report = module.fcatlas_hotspots(12)
    assert "235" in report
    assert "buried" in report.lower()


def test_list_filters(plugin):
    module, _ = plugin
    everything = module.fcatlas_list()
    silenced = module.fcatlas_list("silence")
    assert len(silenced.splitlines()) < len(everything.splitlines())
    assert "LALA" in silenced


def test_compare_offsets_the_second_copy(plugin):
    module, cmd = plugin
    report = module.fcatlas_compare("LALA", "LALA-PG", "1E4K")
    assert len(cmd.of("fetch")) == 2
    assert cmd.of("translate")[0][2]["object"] == "cmp_LALA_PG"
    assert "shared positions" in report
    assert "EU234" in report and "EU235" in report
