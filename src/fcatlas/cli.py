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

    sub.add_parser("validate", help="check every record against real sequence")
    sub.add_parser("structures", help="show structures and their Fc genotypes")

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
        for v in hits:
            print(interface_report(v, cx).summary())
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
