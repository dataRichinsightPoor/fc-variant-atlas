"""Export the atlas to formats other tools can consume.

Everything here is deterministic: no timestamps, no dictionary-order
surprises, no locale-dependent formatting. Running an export twice on the
same input produces byte-identical output, which is what makes the committed
artifacts reviewable in a diff.
"""

from __future__ import annotations

import csv
import io
import json
from typing import Dict, Iterable, List, Optional, Sequence

from .numbering import ISOTYPES, align_isotypes, numbering
from .structure import Complex, classify, load_complexes
from .variants import INTENT_LABELS, Variant, load_variants

CSV_COLUMNS = [
    "id",
    "label",
    "isotype",
    "mutations",
    "disputed_mutations",
    "n_positions",
    "eu_positions",
    "intent",
    "intent_label",
    "imgt",
    "evidence",
    "fcgr_interface_positions",
    "fcrn_interface_positions",
    "off_interface_positions",
    "is_glycan_variant",
    "has_deletion",
    "therapeutics",
    "phenotype",
    "source",
    "notes",
]


def _joined(xs: Iterable) -> str:
    return "; ".join(str(x) for x in xs)


def variant_rows(
    variants: Optional[Sequence[Variant]] = None,
    complexes: Optional[Dict[str, Complex]] = None,
) -> List[dict]:
    """Flat, fully joined rows: curation plus computed structural context."""
    vs = list(variants) if variants is not None else load_variants()
    cx = complexes or load_complexes()
    by_id = {r["id"]: r for r in classify(vs, cx)}
    rows = []
    for v in vs:
        c = by_id[v.id]
        rows.append(
            {
                "id": v.id,
                "label": v.label,
                "isotype": v.isotype,
                "mutations": v.mutation_string(),
                "disputed_mutations": _joined(v.disputed),
                "n_positions": len(v.positions),
                "eu_positions": _joined(v.positions),
                "intent": v.intent,
                "intent_label": INTENT_LABELS.get(v.intent, v.intent),
                "imgt": v.imgt or "",
                "evidence": v.evidence,
                "fcgr_interface_positions": _joined(c["fcgr_interface"]),
                "fcrn_interface_positions": _joined(c["fcrn_interface"]),
                "off_interface_positions": _joined(c["off_interface"]),
                "is_glycan_variant": v.is_glycan_variant,
                "has_deletion": v.has_deletion,
                "therapeutics": _joined(v.therapeutics),
                "phenotype": v.phenotype,
                "source": v.source,
                "notes": " ".join(v.notes.split()),
            }
        )
    return rows


def to_csv(rows: Optional[Sequence[dict]] = None) -> str:
    rs = list(rows) if rows is not None else variant_rows()
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=CSV_COLUMNS, lineterminator="\n")
    w.writeheader()
    for r in rs:
        w.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
    return buf.getvalue()


def to_json(
    variants: Optional[Sequence[Variant]] = None,
    complexes: Optional[Dict[str, Complex]] = None,
) -> str:
    """Full machine-readable atlas, including the alignment and interfaces."""
    vs = list(variants) if variants is not None else load_variants()
    cx = complexes or load_complexes()
    payload = {
        "schema": "fcatlas/1",
        "numbering": "EU (Edelman 1969)",
        "isotypes": {
            iso: {
                "uniprot": _acc(iso),
                "sequence": numbering(iso).sequence,
                "eu_range": [
                    numbering(iso).positions()[0],
                    numbering(iso).positions()[-1],
                ],
                "unnumbered_residue_count": len(numbering(iso).unnumbered_indices),
            }
            for iso in ISOTYPES
        },
        "alignment": [
            {"eu": eu, **col} for eu, col in align_isotypes()
        ],
        "structures": {
            pdb: {
                "partner": c.partner,
                "fc_chains": c.fc_chains,
                "partner_chains": c.partner_chains,
                "fc_glycan_chains": c.fc_glycan_chains,
                "partner_glycan_chains": c.partner_glycan_chains,
                "fc_genotype": c.fc_genotype,
                "wild_type_fc": c.is_wild_type,
                "partner_genotype": c.partner_genotype,
                "wild_type_partner": c.has_wild_type_partner,
                "method": c.method,
                "resolution_angstrom": c.resolution_angstrom,
                "fc_accession": c.fc_accession,
                "fc_accession_database": c.fc_accession_database,
                "fc_reference": c.fc_reference,
                "partner_accessions": c.partner_accessions,
                "partner_accession_databases": c.partner_accession_databases,
                "cutoff_angstrom": c.cutoff,
                "doi": c.doi,
                "caveat": c.caveat,
                "partner_caveat": c.partner_caveat,
                "other_polymers_in_entry": c.other_polymers_in_entry,
                "contacts": {str(k): v for k, v in sorted(c.contacts.items())},
            }
            for pdb, c in sorted(cx.items())
        },
        "variants": variant_rows(vs, cx),
    }
    return json.dumps(payload, indent=1, sort_keys=False, ensure_ascii=True) + "\n"


def _acc(iso: str) -> str:
    from .numbering import ACCESSIONS

    return ACCESSIONS[iso]


def alignment_fasta(isotypes: Sequence[str] = ISOTYPES) -> str:
    """Gapped alignment in FASTA, in EU coordinates, IgG1 as reference."""
    cols = align_isotypes(tuple(isotypes))
    out = []
    for iso in isotypes:
        seq = "".join(col[iso] for _, col in cols)
        out.append(f">{iso}|{_acc(iso)}|EU{cols[0][0]}-{cols[-1][0]}")
        for i in range(0, len(seq), 60):
            out.append(seq[i : i + 60])
    return "\n".join(out) + "\n"


def alignment_block(
    lo: int,
    hi: int,
    variants: Optional[Sequence[Variant]] = None,
    isotypes: Sequence[str] = ISOTYPES,
    width: int = 60,
) -> str:
    """Human-readable aligned block over an EU range, annotated with variants.

    A dot marks identity to IgG1; a dash marks a residue absent from that
    isotype. The annotation line beneath counts how many curated variants
    touch each position, which turns the alignment into a hotspot map.
    """
    vs = list(variants) if variants is not None else load_variants()
    cols = [(eu, col) for eu, col in align_isotypes(tuple(isotypes)) if lo <= eu <= hi]
    if not cols:
        return f"(no aligned positions in EU {lo}-{hi})\n"
    counts: Dict[int, int] = {}
    for v in vs:
        for p in v.positions:
            counts[p] = counts.get(p, 0) + 1
    ref = isotypes[0]
    out: List[str] = []
    for start in range(0, len(cols), width):
        chunk = cols[start : start + width]
        first, last = chunk[0][0], chunk[-1][0]
        out.append(f"EU {first}-{last}")
        for iso in isotypes:
            row = []
            for _, col in chunk:
                aa = col[iso]
                if iso != ref and aa == col[ref] and aa != "-":
                    row.append(".")
                else:
                    row.append(aa)
            out.append(f"  {iso:5s} {''.join(row)}")
        ann = []
        for eu, _ in chunk:
            n = counts.get(eu, 0)
            ann.append(" " if n == 0 else ("+" if n > 9 else str(n)))
        out.append(f"  {'':5s} {''.join(ann)}")
        out.append("")
    return "\n".join(out)


def hotspot_table(variants: Optional[Sequence[Variant]] = None) -> str:
    """Positions ranked by how many curated variants touch them."""
    from .numbering import domain_of
    from .variants import position_index

    vs = list(variants) if variants is not None else load_variants()
    cx = load_complexes()
    idx = position_index(vs)
    rows = sorted(idx.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    n1 = numbering("IgG1")
    out = [
        f"{'EU':>5} {'aa':>3} {'domain':<7} {'n':>3} "
        f"{'closest measured contact':<34} variants"
    ]
    for eu, vlist in rows:
        near = [
            (c.distance(eu), pdb) for pdb, c in cx.items() if c.at_interface(eu)
        ]
        if near:
            d, pdb = min(near)
            receptor = cx[pdb].partner.split(" ")[0].rstrip(",")
            contact = f"{pdb} {d:.2f} A ({receptor})"
        else:
            contact = "-"
        names = ", ".join(v.label for v in vlist[:4])
        if len(vlist) > 4:
            names += f", +{len(vlist) - 4} more"
        out.append(
            f"{eu:>5} {n1.residue(eu) or '-':>3} {domain_of(eu):<7} "
            f"{len(vlist):>3} {contact:<34} {names}"
        )
    return "\n".join(out) + "\n"
