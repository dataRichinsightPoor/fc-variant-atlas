"""Command line interface: python -m fcatlas <command>."""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List, Optional

from . import __version__
from .export import (
    alignment_block,
    alignment_fasta,
    hotspot_table,
    to_csv,
    to_json,
)
from .interface import (
    burial_by_position,
    engineering_against_burial,
    load_interface,
    method_note,
)
from .numbering import ISOTYPES, numbering
from .structure import (
    classify,
    fcrn_provenance,
    interface_report,
    load_complexes,
    pymol_script,
    write_pymol_scripts,
)
from .variants import INTENT_LABELS, find, load_variants, validate_all


def _fmt_variant(v) -> str:
    lines = [
        f"{v.label}  [{v.id}]",
        f"  isotype     {v.isotype}",
        f"  mutations   {v.mutation_string() or '(none: glycan variant)'}",
    ]
    if v.disputed:
        lines.append(
            f"  disputed    {'/'.join(str(s) for s in v.disputed)}"
            "   (see notes)"
        )
    lines += [
        f"  intent      {INTENT_LABELS.get(v.intent, v.intent)}",
        f"  evidence    {v.evidence}",
    ]
    if v.imgt:
        lines.append(f"  IMGT        {v.imgt}")
    if v.aliases:
        lines.append(f"  aliases     {', '.join(v.aliases)}")
    if v.therapeutics:
        lines.append(f"  in clinic   {', '.join(v.therapeutics)}")
    lines.append(f"  phenotype   {v.phenotype}")
    lines.append(f"  source      {v.source}")
    if v.notes:
        lines.append(f"  notes       {' '.join(v.notes.split())}")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="fcatlas",
        description="Aligned, sequence-validated, structure-linked Fc variants.",
    )
    p.add_argument("--version", action="version", version=f"fcatlas {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("list", help="list curated variants")
    s.add_argument("--intent", help="filter by intent")
    s.add_argument("--isotype", choices=ISOTYPES, help="filter by isotype")

    s = sub.add_parser("show", help="show one variant in detail")
    s.add_argument("query", help="id, alias, IMGT code or mutation string")

    s = sub.add_parser("align", help="print the aligned isotype block")
    s.add_argument("--from", dest="lo", type=int, default=231)
    s.add_argument("--to", dest="hi", type=int, default=340)
    s.add_argument("--fasta", action="store_true", help="emit gapped FASTA instead")

    s = sub.add_parser("hotspots", help="positions ranked by variant density")

    s = sub.add_parser("interface", help="test a variant against measured interfaces")
    s.add_argument("query")

    s = sub.add_parser("pymol", help="emit a PyMOL script")
    s.add_argument("query")
    s.add_argument("--pdb", default="1E4K")
    s.add_argument("-o", "--out", help="write here instead of stdout")

    s = sub.add_parser("pymol-all", help="write one PyMOL script per variant")
    s.add_argument("outdir")
    s.add_argument("--pdb", default="1E4K")

    s = sub.add_parser("export", help="write the machine-readable atlas")
    s.add_argument("format", choices=["csv", "json", "fasta"])
    s.add_argument("-o", "--out")

    s = sub.add_parser("burial", help="buried surface area per Fc position")
    s.add_argument("--pdb", default="", help="restrict to one structure")
    s.add_argument("--top", type=int, default=15)
    s.add_argument(
        "--against-density",
        action="store_true",
        help="put engineering frequency beside measured burial",
    )

    s = sub.add_parser("partners", help="named partner residues a variant touches")
    s.add_argument("query")
    s.add_argument("--pdb", default="")

    sub.add_parser("validate", help="check every record against real sequence")
    sub.add_parser("structures", help="show structures and their Fc genotypes")
    sub.add_parser("provenance", help="what is engineered in each deposited entry")

    a = p.parse_args(argv)

    if a.cmd == "list":
        vs = load_variants()
        if a.intent:
            vs = [v for v in vs if v.intent == a.intent]
        if a.isotype:
            vs = [v for v in vs if v.isotype == a.isotype]
        print(f"{'id':22s} {'isotype':7s} {'intent':18s} mutations")
        for v in vs:
            print(
                f"{v.id:22s} {v.isotype:7s} {v.intent:18s} "
                f"{v.mutation_string(mark_disputed=True) or '(glycan variant)'}"
            )
        print(f"\n{len(vs)} variant(s)")
        return 0

    if a.cmd == "show":
        hits = find(a.query)
        if not hits:
            print(f"no variant matching {a.query!r}", file=sys.stderr)
            return 1
        for v in hits:
            print(_fmt_variant(v))
            print()
            print(interface_report(v).summary())
            print()
        return 0

    if a.cmd == "align":
        if a.fasta:
            sys.stdout.write(alignment_fasta())
        else:
            sys.stdout.write(alignment_block(a.lo, a.hi))
        return 0

    if a.cmd == "hotspots":
        sys.stdout.write(hotspot_table())
        return 0

    if a.cmd == "interface":
        hits = find(a.query)
        if not hits:
            print(f"no variant matching {a.query!r}", file=sys.stderr)
            return 1
        cx = load_complexes()
        detail = load_interface()
        for v in hits:
            print(interface_report(v, cx).summary())
            for pdb in sorted(detail):
                rows = [r for r in detail[pdb].for_variant(v) if r.is_buried or r.in_contact]
                if not rows:
                    continue
                print(f"  buried surface in {pdb}:")
                for r in rows:
                    print(f"    {r.summary()}")
                    if r.partners:
                        print(f"        nearest partner residue {r.partners[0]}")
                print(
                    f"    these positions bury {detail[pdb].variant_buried_area(v)} "
                    f"A^2 of the {detail[pdb].total_buried_area} A^2 interface"
                )
            print()
        return 0

    if a.cmd == "pymol":
        hits = find(a.query)
        if not hits:
            print(f"no variant matching {a.query!r}", file=sys.stderr)
            return 1
        v = hits[0]
        stem = f"{v.id}_{a.pdb}"
        text = pymol_script(v, pdb=a.pdb, stem=stem)
        if a.out:
            with open(a.out, "w", encoding="utf-8") as fh:
                fh.write(text)
            print(f"wrote {a.out}")
        else:
            sys.stdout.write(text)
        return 0

    if a.cmd == "pymol-all":
        paths = write_pymol_scripts(a.outdir, pdb=a.pdb)
        print(f"wrote {len(paths)} script(s) to {a.outdir}")
        return 0

    if a.cmd == "export":
        text = {
            "csv": to_csv,
            "json": to_json,
            "fasta": alignment_fasta,
        }[a.format]()
        if a.out:
            os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
            with open(a.out, "w", encoding="utf-8") as fh:
                fh.write(text)
            print(f"wrote {a.out} ({len(text)} bytes)")
        else:
            sys.stdout.write(text)
        return 0

    if a.cmd == "burial":
        detail = load_interface()
        if a.against_density:
            rows = [
                r
                for r in engineering_against_burial()
                if r["records"] or r["delta_sasa"]
            ]
            rows.sort(key=lambda r: (-(r["delta_sasa"] or 0), -r["records"]))
            print("Positions ranked by the largest buried area measured at any")
            print("interface here, with how many curated records change them.")
            print("The two orders are not the same, which is the point of both.")
            header = (
                f"\n{'EU':>5} {'res':>4} {'buried A^2':>11} {'of free':>8} "
                f"{'nearest':>8} {'records':>8}  where"
            )
            print(header)
            for r in rows[: a.top * 2]:
                area = f"{r['delta_sasa']:.1f}" if r["delta_sasa"] else "-"
                frac = (
                    f"{r['buried_fraction'] * 100:.0f}%" if r["buried_fraction"] else "-"
                )
                near = f"{r['min_distance']:.2f}" if r["min_distance"] else "-"
                print(
                    f"{r['eu']:>5} {r['residue'] or '-':>4} {area:>11} {frac:>8} "
                    f"{near:>8} {r['records']:>8}  {r['pdb'] or '-'}"
                )
            return 0
        targets = [a.pdb.upper()] if a.pdb else sorted(detail)
        for pdb in targets:
            if pdb not in detail:
                print(f"no structure {pdb}", file=sys.stderr)
                return 1
            d = detail[pdb]
            print(
                f"{pdb}  {d.total_buried_area} A^2 of Fc surface buried, "
                f"{len(d.buried_positions)} positions, "
                f"glycan {d.glycan_buried_area} A^2"
            )
            for r in d.ranked(a.top):
                print(f"    {r.summary()}")
                # A position can bury 120 A^2 on one heavy chain and 5 on the
                # other. Printing only the representative chain would let a
                # reader treat an asymmetric contact as a symmetric one.
                if r.is_asymmetric:
                    split = ", ".join(
                        f"chain {c} {v['delta_sasa']} A^2"
                        for c, v in sorted(r.per_chain.items())
                    )
                    print(f"        asymmetric between the heavy chains: {split}")
            print()
        note = method_note()
        print(
            f"{note['sasa_algorithm']}, probe {note['probe_radius_angstrom']} A, "
            f"{note['test_points_per_atom']} points per atom, {note['atoms']}"
        )
        print(f"chain choice: {note['chain_choice']}")
        return 0

    if a.cmd == "partners":
        hits = find(a.query)
        if not hits:
            print(f"no variant matching {a.query!r}", file=sys.stderr)
            return 1
        detail = load_interface()
        for v in hits:
            targets = [a.pdb.upper()] if a.pdb else sorted(detail)
            for pdb in targets:
                measured = detail[pdb].for_variant(v)
                rows = [r for r in measured if r.partners]
                if not rows:
                    # An empty table reads exactly like a failed lookup, so say
                    # which of the two silences this is.
                    print(f"{v.label} against {pdb}")
                    if measured:
                        touched = ", ".join(f"EU{r.eu}" for r in measured)
                        print(
                            f"      {touched} lose surface here but no heavy atom"
                            " comes within 5.0 A of the partner"
                        )
                    else:
                        print(
                            "      none of this variant's positions is measured at"
                            f" the interface in {pdb}"
                        )
                    print()
                    continue
                print(f"{v.label} against {pdb}")
                for r in rows:
                    print(f"  EU{r.eu} {r.residue} chain {r.chain}")
                    for partner in r.partners:
                        print(f"      {partner}")
                # T256E is the reason this branch exists: it loses 15.1 A^2 in
                # 4N0U and touches nothing within the cutoff. Dropping it from
                # the output would make the variant look smaller than it is.
                quiet = [r for r in measured if not r.partners]
                for r in quiet:
                    print(f"  EU{r.eu} {r.residue} chain {r.chain}")
                    print(
                        f"      buries {r.delta_sasa} A^2 with no partner atom"
                        " within 5.0 A"
                    )
                print()
        return 0

    if a.cmd == "provenance":
        cx = load_complexes()
        detail = load_interface()
        for pdb, c in sorted(cx.items()):
            print(f"{pdb}  {c.partner}")
            print(
                f"      {c.method}"
                + (f", {c.resolution_angstrom} A" if c.resolution_angstrom else "")
            )
            print(
                f"      Fc {'+'.join(c.fc_chains)} ({c.fc_reference}), "
                f"partner {'+'.join(c.partner_chains)} "
                f"({', '.join(c.partner_accessions) or 'no accession stated'})"
            )
            print(f"      Fc genotype: {c.fc_genotype}")
            print(f"      partner genotype: {c.partner_genotype}")
            for other in c.other_polymers_in_entry:
                print(
                    f"      also in the entry: {other.get('description')} on chain "
                    f"{other.get('chain')}, closest approach to the Fc "
                    f"{other.get('closest_approach_to_fc_angstrom')} A"
                )
            if pdb in detail:
                print(
                    f"      {detail[pdb].total_buried_area} A^2 buried, "
                    f"{len(detail[pdb].contact_positions)} Fc positions in contact"
                )
            if c.caveat:
                print(f"      CAVEAT: {' '.join(c.caveat.split())}")
            if c.partner_caveat:
                print(f"      PARTNER CAVEAT: {' '.join(c.partner_caveat.split())}")
            if c.doi:
                print(f"      {c.doi}")
            print()
        return 0

    if a.cmd == "validate":
        problems = validate_all()
        vs = load_variants()
        n_sub = sum(len(v.substitutions) for v in vs)
        n_dis = sum(len(v.disputed) for v in vs)
        print(f"{len(vs)} records, {n_sub} substitutions validated against sequence")
        print(f"{n_dis} substitution(s) flagged as numbering- or allotype-dependent")
        for iso in ISOTYPES:
            n = numbering(iso)
            print(
                f"  {iso}: {len(n.sequence)} residues, {len(n)} EU-numbered, "
                f"{len(n.unnumbered_indices)} without an EU counterpart"
            )
        if problems:
            print(f"\n{len(problems)} PROBLEM(S):")
            for x in problems:
                print(f"  {x}")
            return 1
        print("\nno problems found")
        return 0

    if a.cmd == "structures":
        cx = load_complexes()
        for pdb, c in sorted(cx.items()):
            print(f"{pdb}  {c.partner}")
            print(f"      Fc chains {'+'.join(c.fc_chains)}, "
                  f"partner {'+'.join(c.partner_chains)}")
            print(f"      Fc genotype: {c.fc_genotype}")
            print(f"      {len(c.contacts)} Fc residues within {c.cutoff} A")
            if c.doi:
                print(f"      {c.doi}")
            if c.caveat:
                print(f"      CAVEAT: {' '.join(c.caveat.split())}")
            print()
        prov = fcrn_provenance()
        print("FcRn structure provenance survey")
        print(f"  {' '.join(prov['_finding'].split())}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
