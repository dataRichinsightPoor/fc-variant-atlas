"""Fc Variant Atlas plugin for PyMOL.

Installs a set of ``fcatlas_*`` commands, and a variant browser if PyMOL was
built with Qt. The commands fetch a deposited complex, colour the measured
interface, mark the positions a chosen variant changes, and print what the
measurement rests on before drawing anything: the resolution, the genotype of
the Fc chain in that entry, the genotype of the partner, and how much surface
each position actually gives up when the partner binds.

Three deliberate choices, each of which is a claim about the picture:

Sugars are drawn as spheres and never as lines. Glycan bond connectivity in a
deposited entry is inferred by the viewer, not deposited, and a line drawing
asserts bonds that are not in the coordinate file.

Burial is loaded into the B-factor column, so any PyMOL colouring or selection
that reads B reads measured buried area. Nothing else is written there.

Every scene prints the genotype of the structure it just fetched. The Fc in
the FcRn complex carries YTE, and the receptor in the only FcgammaRIIIa entry
carries four substitutions, so a scene that stays silent about provenance
invites a reader to mistake an engineered molecule for a reference one.

Install: PyMOL, Plugin, Plugin Manager, Install New Plugin, and choose the
zip file, or the fcatlas_pymol directory. Then type ``fcatlas`` in the PyMOL
prompt for the command list.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence

__version__ = "1.1.0"
__author__ = "Ermelinda Damko"

REPO_URL = "https://github.com/dataRichinsightPoor/fc-variant-atlas"

FC_COLOR = "grey90"
PARTNER_COLOR = "grey70"
INTERFACE_COLOR = "palecyan"
VARIANT_COLOR = "red"
GLYCAN_COLOR = "gold"
BURIAL_PALETTE = "white_yellow_orange_red"

_CACHE: Dict[str, dict] = {}


# ----------------------------------------------------------------------------
# data loading
# ----------------------------------------------------------------------------


def _candidate_dirs() -> List[Path]:
    """Where to look for atlas.json and interface_detail.json, in order."""
    found: List[Path] = []
    env = os.environ.get("FCATLAS_DATA")
    if env:
        found.append(Path(env))
    try:  # an installed fcatlas package is the best source
        from importlib import resources

        found.append(Path(str(resources.files("fcatlas").joinpath("data"))))
    except Exception:
        pass
    here = Path(__file__).resolve().parent
    found.append(here / "data")
    for parent in list(here.parents)[:5]:
        found.append(parent / "data")
    return found


def data_dir() -> Path:
    """First directory that holds the atlas export."""
    for candidate in _candidate_dirs():
        if (candidate / "atlas.json").is_file():
            return candidate
    searched = "\n  ".join(str(c) for c in _candidate_dirs())
    raise FileNotFoundError(
        "Could not find atlas.json. Set FCATLAS_DATA to the atlas data "
        f"directory, or install the fcatlas package. Looked in:\n  {searched}"
    )


def _load(name: str) -> dict:
    if name not in _CACHE:
        path = data_dir() / name
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found in {path.parent}")
        _CACHE[name] = json.loads(path.read_text(encoding="utf-8"))
    return _CACHE[name]


def atlas() -> dict:
    return _load("atlas.json")


def interface_detail() -> dict:
    """Per-position areas and partner residues; optional."""
    try:
        return _load("interface_detail.json")
    except FileNotFoundError:
        return {"complexes": {}}


def variants() -> List[dict]:
    return atlas()["variants"]


def structures() -> Dict[str, dict]:
    return atlas()["structures"]


def find_variant(variant_id: str) -> dict:
    """Look a variant up by id, label or alias, case-insensitively."""
    wanted = str(variant_id).strip().lower()
    for record in variants():
        if wanted in {
            str(record.get("id", "")).lower(),
            str(record.get("label", "")).lower(),
        }:
            return record
    loose = [r for r in variants() if wanted in str(r.get("id", "")).lower()]
    if len(loose) == 1:
        return loose[0]
    if loose:
        names = ", ".join(r["id"] for r in loose[:12])
        raise KeyError(f"{variant_id!r} matches several variants: {names}")
    raise KeyError(
        f"No variant {variant_id!r}. Run fcatlas_list to see the {len(variants())} ids."
    )


def positions_of(record: dict) -> List[int]:
    field = record.get("eu_positions", "")
    return [int(x) for x in str(field).replace(",", ";").split(";") if x.strip()]


def default_structure(record: dict) -> str:
    """FcRn complex for half-life variants, otherwise the FcgammaR complex."""
    if record.get("fcrn_interface_positions") and not record.get(
        "fcgr_interface_positions"
    ):
        return "4N0U"
    return "1E4K"


# ----------------------------------------------------------------------------
# reporting, all pure text so it can be tested without PyMOL
# ----------------------------------------------------------------------------


def provenance_lines(pdb: str) -> List[str]:
    """What the structure is, and what is engineered in it."""
    record = structures().get(pdb)
    if record is None:
        return [f"No structure {pdb} in the atlas."]
    lines = [
        f"{pdb}: Fc plus {record.get('partner', 'partner')}",
        f"   {record.get('method', 'method not stated')}"
        + (
            f", {record['resolution_angstrom']} A"
            if record.get("resolution_angstrom")
            else ""
        ),
        f"   Fc chains {'+'.join(record.get('fc_chains', []))}"
        + (
            f", Fc glycan chains {'+'.join(record['fc_glycan_chains'])}"
            if record.get("fc_glycan_chains")
            else ""
        )
        + f", partner chains {'+'.join(record.get('partner_chains', []))}",
        f"   Fc genotype: {record.get('fc_genotype', 'not recorded')}",
        f"   partner genotype: {record.get('partner_genotype', 'not recorded')}",
    ]
    if record.get("caveat"):
        lines.append(f"   caveat: {record['caveat']}")
    if record.get("partner_caveat"):
        lines.append(f"   partner caveat: {record['partner_caveat']}")
    for other in record.get("other_polymers_in_entry", []):
        distance = other.get("closest_approach_to_fc_angstrom")
        lines.append(
            f"   also in this entry: {other.get('description')} on chain "
            f"{other.get('chain')}"
            + (f", closest approach to the Fc {distance} A" if distance else "")
        )
    if record.get("doi"):
        lines.append(f"   source: {record['doi']}")
    return lines


def position_rows(pdb: str, eus: Optional[Sequence[int]] = None) -> List[dict]:
    """Measured detail for the requested positions, most buried first."""
    complexes = interface_detail().get("complexes", {})
    detail = complexes.get(pdb, {}).get("positions", {})
    contacts = structures().get(pdb, {}).get("contacts", {})
    wanted = list(eus) if eus is not None else sorted(
        {int(k) for k in detail} | {int(k) for k in contacts}
    )
    rows = []
    for eu in wanted:
        record = detail.get(str(eu), {})
        # The distance reported belongs to the same heavy chain as the area, so
        # the row describes one copy of the residue. Only fall back to the
        # cross-chain contact table where there is no measurement at all.
        distance = record.get("min_distance")
        if distance is None and not record:
            distance = contacts.get(str(eu))
        rows.append(
            {
                "eu": int(eu),
                "residue": record.get("residue"),
                "chain": record.get("chain"),
                "min_distance": distance,
                "delta_sasa": record.get("delta_sasa"),
                "buried_fraction": record.get("buried_fraction"),
                "partners": record.get("partners", []),
                "per_chain": record.get("per_chain", {}),
                "measured": bool(record),
            }
        )
    rows.sort(key=lambda r: -(r["delta_sasa"] or 0))
    return rows


def interface_lines(variant_id: str, pdb: Optional[str] = None) -> List[str]:
    """A table of what this variant's positions do at this interface."""
    record = find_variant(variant_id)
    pdb = pdb or default_structure(record)
    eus = positions_of(record)
    rows = position_rows(pdb, eus)
    lines = [
        f"{record['id']} ({record.get('isotype', '')}) {record.get('mutations', '')}"
        f" against {pdb}",
        "   EU  res  chain  nearest  buried A^2   of free  closest partner residue",
    ]
    total = 0.0
    for row in rows:
        area = row["delta_sasa"]
        total += area or 0.0
        nearest = f"{row['min_distance']:.2f}" if row["min_distance"] else "  -  "
        fraction = f"{row['buried_fraction'] * 100:3.0f}%" if row["buried_fraction"] else "  - "
        partner = row["partners"][0] if row["partners"] else None
        partner_text = (
            f"{partner['residue']}{partner['seq']} chain {partner['chain']}"
            + (" (glycan)" if partner.get("is_glycan") else "")
            if partner
            else "no contact within the cutoff"
        )
        lines.append(
            f"   {row['eu']:>3}  {row['residue'] or '?':>3}  {row['chain'] or '?':>4}   "
            f"{nearest:>6}  {(area if area is not None else 0):>9.1f}   {fraction:>6}  "
            f"{partner_text}"
        )
        asymmetry = row["per_chain"]
        if len(asymmetry) > 1:
            values = {c: v.get("delta_sasa", 0.0) for c, v in asymmetry.items()}
            low, high = min(values.values()), max(values.values())
            if high >= 10 and low <= 0.25 * high:
                spread = ", ".join(f"chain {c} {v:.1f}" for c, v in values.items())
                lines.append(
                    f"        this position is asymmetric between the heavy chains: {spread}"
                )
    lines.append(f"   total surface buried by these positions: {total:.1f} A^2")
    unmeasured = [r["eu"] for r in rows if not r["measured"]]
    if unmeasured:
        lines.append(
            "   no measured burial at "
            + ", ".join(f"EU{eu}" for eu in unmeasured)
            + f": these positions are not at the interface in {pdb}"
        )
    return lines


def variant_lines(record: dict) -> List[str]:
    """The curated record, without the structure."""
    lines = [
        f"{record['id']}  {record.get('isotype', '')}  {record.get('mutations', '')}",
        f"   intent: {record.get('intent_label', record.get('intent', ''))}",
    ]
    if record.get("imgt"):
        lines.append(f"   IMGT designation: {record['imgt']}")
    if record.get("disputed_mutations"):
        lines.append(
            f"   numbering dependent: also written {record['disputed_mutations']} "
            "under the other convention"
        )
    if record.get("phenotype"):
        lines.append(f"   phenotype: {record['phenotype']}")
    if record.get("therapeutics"):
        lines.append(f"   in: {record['therapeutics']}")
    if record.get("source"):
        lines.append(f"   source: {record['source']}")
    return lines


def list_lines(intent: str = "", isotype: str = "") -> List[str]:
    rows = variants()
    if intent:
        rows = [r for r in rows if str(r.get("intent", "")).lower() == intent.lower()]
    if isotype:
        rows = [r for r in rows if str(r.get("isotype", "")).lower() == isotype.lower()]
    lines = [f"{len(rows)} variants" + (f" [intent={intent}]" if intent else "")]
    for record in rows:
        lines.append(
            f"   {record['id']:<22} {record.get('isotype', ''):<5} "
            f"{record.get('mutations', '')}"
        )
    return lines


def hotspot_lines(limit: int = 20) -> List[str]:
    """Positions ranked by how often they are engineered, with burial beside."""
    counts: Dict[int, int] = {}
    for record in variants():
        for eu in set(positions_of(record)):
            counts[eu] = counts.get(eu, 0) + 1
    best: Dict[int, tuple] = {}
    for pdb, entry in interface_detail().get("complexes", {}).items():
        for eu, record in entry.get("positions", {}).items():
            area = record.get("delta_sasa") or 0.0
            if int(eu) not in best or area > best[int(eu)][0]:
                best[int(eu)] = (area, pdb)
    lines = [
        "positions ranked by how many curated records change them, with the largest",
        "buried area measured at any interface in this set beside them",
        "   EU  records  buried A^2  where",
    ]
    for eu, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]:
        area, pdb = best.get(eu, (None, "not at any interface here"))
        area_text = f"{area:9.1f}" if area else "        -"
        lines.append(f"   {eu:>3}  {count:>7}  {area_text}  {pdb}")
    return lines


# ----------------------------------------------------------------------------
# PyMOL commands
# ----------------------------------------------------------------------------


def _cmd():
    from pymol import cmd

    return cmd


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(name))


def _report(lines: Sequence[str]) -> str:
    text = "\n".join(lines)
    print(text)
    return text


def fcatlas(*_args) -> str:
    """Print the command list."""
    return _report(
        [
            f"Fc Variant Atlas {__version__} — {len(variants())} curated variants",
            f"data: {data_dir()}",
            f"repository: {REPO_URL}",
            "",
            "   fcatlas_list [intent [, isotype]]      list variants",
            "   fcatlas_show ID [, pdb]                fetch and build a scene",
            "   fcatlas_interface ID [, pdb]           measured contacts and burial",
            "   fcatlas_partners ID [, pdb]            show the partner residues too",
            "   fcatlas_burial [pdb]                   colour by buried area",
            "   fcatlas_compare ID1, ID2 [, pdb]       two variants side by side",
            "   fcatlas_provenance [pdb]               what is engineered in the entry",
            "   fcatlas_hotspots [limit]               engineering against burial",
            "",
            "   intents: " + ", ".join(sorted({r.get("intent", "") for r in variants()})),
            "   structures: " + ", ".join(structures()),
        ]
    )


def fcatlas_list(intent: str = "", isotype: str = "") -> str:
    """List curated variants, optionally filtered by intent and isotype."""
    return _report(list_lines(intent, isotype))


def fcatlas_provenance(pdb: str = "") -> str:
    """Print what each structure is and what is engineered in it."""
    targets = [pdb] if pdb else list(structures())
    lines: List[str] = []
    for target in targets:
        lines.extend(provenance_lines(target))
    return _report(lines)


def fcatlas_interface(variant_id: str, pdb: str = "") -> str:
    """Print measured contacts and buried area for one variant."""
    return _report(interface_lines(variant_id, pdb or None))


def fcatlas_hotspots(limit: int = 20) -> str:
    """Print positional density beside measured burial."""
    return _report(hotspot_lines(int(limit)))


def fcatlas_show(
    variant_id: str,
    pdb: str = "",
    fetch: int = 1,
    keep: int = 0,
    surface: int = 0,
) -> str:
    """Fetch a complex and build a scene for one variant.

    variant_id  a curated id, for example LALA-PG
    pdb         which complex, default 1E4K for effector variants and 4N0U
                for half-life variants
    fetch       0 to reuse an object already loaded under the structure name
    keep        1 to leave existing objects in the session
    surface     1 to draw the partner as a surface instead of cartoon
    """
    cmd = _cmd()
    record = find_variant(variant_id)
    pdb = (pdb or default_structure(record)).upper()
    structure = structures().get(pdb)
    if structure is None:
        raise KeyError(f"No structure {pdb}. Known: {', '.join(structures())}")

    name = f"fc_{_safe(record['id'])}_{pdb}"
    eus = positions_of(record)
    fc_sel = "+".join(structure["fc_chains"])
    partner_sel = "+".join(structure["partner_chains"])
    glycan_chains = structure.get("fc_glycan_chains", []) + structure.get(
        "partner_glycan_chains", []
    )

    if not int(keep):
        cmd.delete("all")
    if int(fetch):
        cmd.fetch(pdb, name, type="cif", async_=0)
    elif name not in cmd.get_names("objects"):
        cmd.set_name(pdb, name)

    cmd.remove(f"{name} and solvent")
    cmd.hide("everything", name)
    cmd.show("cartoon", f"{name} and chain {fc_sel} and polymer")
    if int(surface):
        cmd.show("surface", f"{name} and chain {partner_sel} and polymer")
        cmd.set("transparency", 0.25, name)
    else:
        cmd.show("cartoon", f"{name} and chain {partner_sel} and polymer")
    cmd.color(FC_COLOR, f"{name} and chain {fc_sel}")
    cmd.color(PARTNER_COLOR, f"{name} and chain {partner_sel}")

    interface_eus = sorted(int(k) for k in structure.get("contacts", {}))
    if interface_eus:
        selection = "+".join(str(x) for x in interface_eus)
        cmd.select(
            f"{name}_interface", f"{name} and chain {fc_sel} and resi {selection}"
        )
        cmd.color(INTERFACE_COLOR, f"{name}_interface")
        cmd.show("sticks", f"{name}_interface and sidechain")

    if eus:
        selection = "+".join(str(x) for x in eus)
        cmd.select(f"{name}_variant", f"{name} and chain {fc_sel} and resi {selection}")
        cmd.color(VARIANT_COLOR, f"{name}_variant")
        cmd.show("sticks", f"{name}_variant")
        cmd.set("stick_radius", 0.22, f"{name}_variant")
        cmd.label(f"{name}_variant and name CA", "'%s%s' % (resn, resi)")

    # Sugars as spheres: sphere rendering makes no bond-connectivity claim.
    if glycan_chains:
        cmd.select(
            f"{name}_glycans",
            f"{name} and chain {'+'.join(glycan_chains)} and not polymer",
        )
    else:
        cmd.select(f"{name}_glycans", f"{name} and organic and not polymer")
    if cmd.count_atoms(f"{name}_glycans"):
        cmd.show("spheres", f"{name}_glycans")
        cmd.set("sphere_scale", 0.32, f"{name}_glycans")
        cmd.color(GLYCAN_COLOR, f"{name}_glycans")

    cmd.deselect()
    cmd.orient(name)
    if eus and cmd.count_atoms(f"{name}_variant"):
        cmd.zoom(f"{name}_variant", 8)

    lines = variant_lines(record) + [""] + provenance_lines(pdb) + [""]
    lines += interface_lines(record["id"], pdb)
    lines += [
        "",
        f"objects: {name}, selections {name}_interface, {name}_variant, {name}_glycans",
        "sugars are drawn as spheres because their bond connectivity is inferred "
        "rather than deposited",
        "next: fcatlas_burial to colour by measured buried area, or fcatlas_partners "
        "to show the partner side",
    ]
    return _report(lines)


def fcatlas_burial(pdb: str = "", obj: str = "") -> str:
    """Write measured buried area into B-factors and colour by it.

    Nothing else is written to the B-factor column, so any selection on b is a
    selection on measured buried surface area.
    """
    cmd = _cmd()
    names = cmd.get_names("objects")
    if not names:
        return _report(["Nothing loaded. Run fcatlas_show first."])
    target = obj or names[0]
    if not pdb:
        pdb = next((p for p in structures() if p in target.upper()), "")
    if not pdb:
        return _report(
            [f"Could not tell which structure {target} is. Pass the pdb id."]
        )
    pdb = pdb.upper()

    structure = structures()[pdb]
    fc_sel = "+".join(structure["fc_chains"])
    rows = [r for r in position_rows(pdb) if r["delta_sasa"]]
    cmd.alter(target, "b=0.0")
    for row in rows:
        chain = row["chain"] or structure["fc_chains"][0]
        cmd.alter(
            f"{target} and chain {chain} and resi {row['eu']}",
            f"b={row['delta_sasa']}",
        )
    cmd.sort(target)
    ceiling = max((r["delta_sasa"] or 0) for r in rows) if rows else 1.0
    cmd.spectrum(
        "b", BURIAL_PALETTE, f"{target} and chain {fc_sel}", minimum=0, maximum=ceiling
    )
    top = rows[:8]
    lines = [
        f"buried area written to B-factors for {target} ({pdb}), "
        f"{len(rows)} positions, scale 0 to {ceiling:.0f} A^2",
        "   select deeply_buried, b > 50   to work with the positions that give up "
        "the most surface",
        "   most buried here:",
    ]
    for row in top:
        lines.append(
            f"      EU{row['eu']} {row['residue']} chain {row['chain']}: "
            f"{row['delta_sasa']:.1f} A^2, {row['buried_fraction'] * 100:.0f}% of its "
            "free surface"
        )
    return _report(lines)


def fcatlas_partners(variant_id: str, pdb: str = "", cutoff: float = 5.0) -> str:
    """Show the partner residues that contact this variant's positions."""
    cmd = _cmd()
    record = find_variant(variant_id)
    pdb = (pdb or default_structure(record)).upper()
    structure = structures()[pdb]
    name = f"fc_{_safe(record['id'])}_{pdb}"
    if name not in cmd.get_names("objects"):
        fcatlas_show(record["id"], pdb)

    rows = position_rows(pdb, positions_of(record))
    partner_residues = sorted(
        {(p["chain"], p["seq"]) for row in rows for p in row["partners"]}
    )
    if not partner_residues:
        return _report(
            [
                f"{record['id']} has no position in contact with the partner in {pdb} "
                f"within {cutoff} A."
            ]
        )

    by_chain: Dict[str, List[int]] = {}
    for chain, seq in partner_residues:
        by_chain.setdefault(chain, []).append(seq)
    clauses = [
        f"(chain {chain} and resi {'+'.join(str(s) for s in sorted(seqs))})"
        for chain, seqs in sorted(by_chain.items())
    ]
    cmd.select(f"{name}_partner_side", f"{name} and (" + " or ".join(clauses) + ")")
    cmd.show("sticks", f"{name}_partner_side and sidechain")
    cmd.color("skyblue", f"{name}_partner_side")
    cmd.deselect()

    lines = [
        f"{record['id']} against {pdb}: partner residues within {cutoff} A of the "
        "changed positions",
    ]
    for row in rows:
        if not row["partners"]:
            continue
        listed = ", ".join(
            f"{p['residue']}{p['seq']} chain {p['chain']} {p['distance']:.2f} A"
            + (" (glycan)" if p.get("is_glycan") else "")
            for p in row["partners"][:6]
        )
        lines.append(f"   EU{row['eu']} {row['residue']}: {listed}")
    lines.append(f"selection {name}_partner_side holds the partner side")
    return _report(lines)


def fcatlas_compare(
    variant_a: str, variant_b: str, pdb: str = "", offset: float = 45.0
) -> str:
    """Load one complex twice and mark a different variant on each copy."""
    cmd = _cmd()
    first, second = find_variant(variant_a), find_variant(variant_b)
    pdb = (pdb or default_structure(first)).upper()
    structure = structures()[pdb]
    fc_sel = "+".join(structure["fc_chains"])

    cmd.delete("all")
    names = []
    for index, record in enumerate((first, second)):
        name = f"cmp_{_safe(record['id'])}"
        names.append(name)
        cmd.fetch(pdb, name, type="cif", async_=0)
        cmd.remove(f"{name} and solvent")
        cmd.hide("everything", name)
        cmd.show("cartoon", f"{name} and polymer")
        cmd.color(FC_COLOR, name)
        cmd.color(PARTNER_COLOR, f"{name} and chain {'+'.join(structure['partner_chains'])}")
        eus = positions_of(record)
        if eus:
            selection = "+".join(str(x) for x in eus)
            cmd.select(
                f"{name}_variant", f"{name} and chain {fc_sel} and resi {selection}"
            )
            cmd.color(VARIANT_COLOR if index == 0 else "marine", f"{name}_variant")
            cmd.show("sticks", f"{name}_variant")
        if index == 1:
            cmd.translate([float(offset), 0.0, 0.0], object=name)
    cmd.deselect()
    cmd.zoom("all", 5)

    shared = sorted(set(positions_of(first)) & set(positions_of(second)))
    rows_a = position_rows(pdb, positions_of(first))
    rows_b = position_rows(pdb, positions_of(second))
    area_a = sum(r["delta_sasa"] or 0 for r in rows_a)
    area_b = sum(r["delta_sasa"] or 0 for r in rows_b)
    lines = [
        f"{first['id']} in red on the left, {second['id']} in blue on the right, "
        f"both on {pdb}",
        f"   {first['id']}: {first.get('mutations', '')}, "
        f"{area_a:.1f} A^2 of interface touched",
        f"   {second['id']}: {second.get('mutations', '')}, "
        f"{area_b:.1f} A^2 of interface touched",
        "   shared positions: "
        + (", ".join(f"EU{eu}" for eu in shared) if shared else "none"),
    ]
    lines.extend(provenance_lines(pdb))
    return _report(lines)


COMMANDS = {
    "fcatlas": fcatlas,
    "fcatlas_list": fcatlas_list,
    "fcatlas_show": fcatlas_show,
    "fcatlas_interface": fcatlas_interface,
    "fcatlas_partners": fcatlas_partners,
    "fcatlas_burial": fcatlas_burial,
    "fcatlas_compare": fcatlas_compare,
    "fcatlas_provenance": fcatlas_provenance,
    "fcatlas_hotspots": fcatlas_hotspots,
}


def register_commands() -> List[str]:
    """Make the fcatlas_* commands available at the PyMOL prompt."""
    from pymol import cmd

    for name, function in COMMANDS.items():
        cmd.extend(name, function)
    try:
        cmd.auto_arg[0].update(
            {
                name: [lambda: [r["id"] for r in variants()], "variant", ", "]
                for name in ("fcatlas_show", "fcatlas_interface", "fcatlas_partners")
            }
        )
    except Exception:
        pass
    return sorted(COMMANDS)


# ----------------------------------------------------------------------------
# Qt browser
# ----------------------------------------------------------------------------


def make_dialog(parent=None):
    """Build the variant browser. Requires a Qt-enabled PyMOL."""
    from pymol.Qt import QtWidgets

    dialog = QtWidgets.QDialog(parent)
    dialog.setWindowTitle(f"Fc Variant Atlas {__version__}")
    dialog.resize(760, 560)
    layout = QtWidgets.QVBoxLayout(dialog)

    filters = QtWidgets.QHBoxLayout()
    search = QtWidgets.QLineEdit()
    search.setPlaceholderText("filter by id, mutation or phenotype")
    intent_box = QtWidgets.QComboBox()
    intent_box.addItem("every intent", "")
    for intent in sorted({r.get("intent", "") for r in variants()}):
        intent_box.addItem(intent, intent)
    structure_box = QtWidgets.QComboBox()
    for pdb, record in structures().items():
        structure_box.addItem(f"{pdb} — {record.get('partner', '')}", pdb)
    filters.addWidget(search, 3)
    filters.addWidget(intent_box, 1)
    filters.addWidget(structure_box, 2)
    layout.addLayout(filters)

    table = QtWidgets.QTableWidget(0, 4)
    table.setHorizontalHeaderLabels(["variant", "isotype", "substitutions", "intent"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    layout.addWidget(table, 3)

    report = QtWidgets.QPlainTextEdit()
    report.setReadOnly(True)
    layout.addWidget(report, 2)

    buttons = QtWidgets.QHBoxLayout()
    show_button = QtWidgets.QPushButton("Show scene")
    burial_button = QtWidgets.QPushButton("Colour by buried area")
    partners_button = QtWidgets.QPushButton("Show partner side")
    provenance_button = QtWidgets.QPushButton("Provenance")
    for widget in (show_button, burial_button, partners_button, provenance_button):
        buttons.addWidget(widget)
    layout.addLayout(buttons)

    rows: List[dict] = []

    def refresh():
        rows.clear()
        needle = search.text().strip().lower()
        intent = intent_box.currentData()
        for record in variants():
            if intent and record.get("intent") != intent:
                continue
            haystack = " ".join(
                str(record.get(key, ""))
                for key in ("id", "label", "mutations", "phenotype", "therapeutics", "imgt")
            ).lower()
            if needle and needle not in haystack:
                continue
            rows.append(record)
        table.setRowCount(len(rows))
        for index, record in enumerate(rows):
            for column, key in enumerate(("id", "isotype", "mutations", "intent_label")):
                table.setItem(
                    index, column, QtWidgets.QTableWidgetItem(str(record.get(key, "")))
                )
        table.resizeColumnsToContents()

    def selected() -> Optional[dict]:
        index = table.currentRow()
        if 0 <= index < len(rows):
            return rows[index]
        return None

    def describe():
        record = selected()
        if record:
            report.setPlainText(
                "\n".join(
                    variant_lines(record)
                    + [""]
                    + interface_lines(record["id"], structure_box.currentData())
                )
            )

    def run(function, *args):
        record = selected()
        if not record:
            report.setPlainText("Select a variant first.")
            return
        try:
            report.setPlainText(function(record["id"], *args))
        except Exception as error:  # keep the dialog alive on a bad scene
            report.setPlainText(f"{type(error).__name__}: {error}")

    search.textChanged.connect(refresh)
    intent_box.currentIndexChanged.connect(refresh)
    table.itemSelectionChanged.connect(describe)
    structure_box.currentIndexChanged.connect(describe)
    show_button.clicked.connect(lambda: run(fcatlas_show, structure_box.currentData()))
    partners_button.clicked.connect(
        lambda: run(fcatlas_partners, structure_box.currentData())
    )
    burial_button.clicked.connect(
        lambda: report.setPlainText(fcatlas_burial(structure_box.currentData()))
    )
    provenance_button.clicked.connect(
        lambda: report.setPlainText(fcatlas_provenance(structure_box.currentData()))
    )

    refresh()
    return dialog


_DIALOG = None


def run_plugin_gui():
    """Open the browser, or fall back to the command list without Qt."""
    global _DIALOG
    try:
        if _DIALOG is None:
            _DIALOG = make_dialog()
        _DIALOG.show()
        _DIALOG.raise_()
        return _DIALOG
    except Exception as error:
        print(
            f"The Fc Variant Atlas browser needs a Qt-enabled PyMOL ({error}). "
            "The fcatlas_* commands work regardless; type fcatlas for the list."
        )
        return fcatlas()


def __init_plugin__(app=None):
    """PyMOL plugin entry point."""
    registered = register_commands()
    try:
        from pymol.plugins import addmenuitemqt

        addmenuitemqt("Fc Variant Atlas", run_plugin_gui)
    except Exception:
        pass
    print(
        f"Fc Variant Atlas {__version__} loaded: {len(registered)} commands, "
        f"{len(variants())} variants. Type fcatlas for the list."
    )


try:  # so `run fcatlas_pymol/__init__.py` also works
    register_commands()
except Exception:
    pass
