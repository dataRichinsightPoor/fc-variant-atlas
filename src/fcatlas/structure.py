"""Structural context for Fc variants.

Two things happen here. First, each variant position is tested against
measured interface contacts from deposited co-crystal structures, so the
question "is this substitution actually where the receptor binds?" has a
numeric answer instead of a narrative one. Second, ready-to-run PyMOL scripts
and browser viewer scenes are generated, so a variant can be looked at
without any manual selection-string writing.

The contact distances are minimum heavy-atom distances from the deposited
coordinates, computed once and committed to data/structures.json rather than
recomputed at import, which keeps the package importable without network
access and keeps the numbers reproducible.

A caveat is carried with the data rather than left in a paper: no deposited
human Fc-FcRn co-crystal structure contains a wild-type IgG1 Fc. The
structure used here for the FcRn interface, 4N0U, carries YTE, and the YTE
substitutions lie inside the contact patch being measured.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib import resources
from typing import Dict, Iterable, List, Optional, Sequence

from .variants import Variant, load_variants

#: Default per-structure colours used by the generated scenes.
PARTNER_COLOR = "grey70"
FC_COLOR = "grey90"
INTERFACE_COLOR = "palecyan"
VARIANT_COLOR = "red"


@dataclass(frozen=True)
class Complex:
    """One deposited structure with a measured Fc interface."""

    pdb: str
    partner: str
    partner_chains: List[str]
    fc_chains: List[str]
    fc_genotype: str
    cutoff: float
    contacts: Dict[int, float]
    doi: str = ""
    caveat: str = ""
    partner_genotype: str = ""
    partner_caveat: str = ""
    method: str = ""
    resolution_angstrom: Optional[float] = None
    fc_accession: str = ""
    fc_accession_database: str = ""
    partner_accessions: List[str] = field(default_factory=list)
    partner_accession_databases: List[str] = field(default_factory=list)
    fc_glycan_chains: List[str] = field(default_factory=list)
    partner_glycan_chains: List[str] = field(default_factory=list)
    deposited_fc_mutation: Optional[str] = None
    deposited_partner_mutation: Optional[str] = None
    other_polymers_in_entry: List[dict] = field(default_factory=list)

    @property
    def is_wild_type(self) -> bool:
        return self.fc_genotype.strip().lower() == "wild type"

    @property
    def has_wild_type_partner(self) -> bool:
        """Whether the binding partner in this entry is unengineered.

        The provenance problem is symmetric. An entry can carry a wild-type Fc
        and still be an engineered picture of the interface, because the
        receptor was mutated to make the complex crystallize.
        """
        return self.deposited_partner_mutation is None

    @property
    def interface_positions(self) -> List[int]:
        return sorted(self.contacts)

    def distance(self, eu: int) -> Optional[float]:
        """Minimum heavy-atom distance from EU position to the partner."""
        return self.contacts.get(eu)

    def at_interface(self, eu: int) -> bool:
        return eu in self.contacts

    @property
    def fc_reference(self) -> str:
        """Accession with the database that issued it."""
        if not self.fc_accession:
            return ""
        return f"{self.fc_accession_database or 'unstated'}:{self.fc_accession}"


def _payload() -> dict:
    return json.loads(
        resources.files("fcatlas")
        .joinpath("data", "structures.json")
        .read_text(encoding="utf-8")
    )


def load_complexes() -> Dict[str, Complex]:
    """Load the committed interface measurements."""
    data = _payload()["complexes"]
    out: Dict[str, Complex] = {}
    for pdb, rec in data.items():
        out[pdb] = Complex(
            pdb=pdb,
            partner=rec["partner"],
            partner_chains=list(rec["partner_chains"]),
            fc_chains=list(rec["fc_chains"]),
            fc_genotype=rec["fc_genotype"],
            cutoff=float(rec["cutoff"]),
            contacts={int(k): float(v) for k, v in rec["contacts"].items()},
            doi=rec.get("doi", ""),
            caveat=rec.get("caveat", ""),
            partner_genotype=rec.get("partner_genotype", ""),
            partner_caveat=rec.get("partner_caveat", ""),
            method=rec.get("method", ""),
            resolution_angstrom=rec.get("resolution_angstrom"),
            fc_accession=rec.get("fc_accession", ""),
            fc_accession_database=rec.get("fc_accession_database") or "",
            partner_accessions=list(rec.get("partner_accessions", [])),
            partner_accession_databases=[
                d or "" for d in rec.get("partner_accession_databases", [])
            ],
            fc_glycan_chains=list(rec.get("fc_glycan_chains", [])),
            partner_glycan_chains=list(rec.get("partner_glycan_chains", [])),
            deposited_fc_mutation=rec.get("deposited_fc_mutation"),
            deposited_partner_mutation=rec.get("deposited_partner_mutation"),
            other_polymers_in_entry=list(rec.get("other_polymers_in_entry", [])),
        )
    return out


def fcrn_provenance() -> dict:
    """The FcRn structure-provenance survey and its finding."""
    return _payload()["fcrn_provenance_survey"]


# ------------------------------------------------------------------ interface Q
@dataclass
class InterfaceReport:
    """Where a variant's positions sit relative to measured interfaces."""

    variant: str
    rows: List[dict]

    def touching(self, pdb: str) -> List[int]:
        return [r["eu"] for r in self.rows if r["distances"].get(pdb) is not None]

    def summary(self) -> str:
        lines = [f"{self.variant}: {len(self.rows)} position(s)"]
        for r in self.rows:
            hits = [
                f"{pdb} {d:.2f} A"
                for pdb, d in sorted(r["distances"].items())
                if d is not None
            ]
            where = "; ".join(hits) if hits else "no measured interface contact"
            lines.append(f"  EU {r['eu']} ({r['change']}): {where}")
        return "\n".join(lines)


def interface_report(
    variant: Variant, complexes: Optional[Dict[str, Complex]] = None
) -> InterfaceReport:
    """Test every position of a variant against every measured interface."""
    cx = complexes or load_complexes()
    rows = []
    changes = {s.eu: str(s) for s in variant.substitutions}
    changes.update({d.eu: str(d) for d in variant.deletions})
    changes.update({s.eu: f"{s} (disputed)" for s in variant.disputed})
    for eu in variant.positions:
        rows.append(
            {
                "eu": eu,
                "change": changes.get(eu, ""),
                "distances": {pdb: c.distance(eu) for pdb, c in cx.items()},
            }
        )
    return InterfaceReport(variant.label, rows)


def classify(
    variants: Optional[Sequence[Variant]] = None,
    complexes: Optional[Dict[str, Complex]] = None,
) -> List[dict]:
    """Summarize, per variant, which measured interfaces its positions touch.

    This is the table that separates variants acting inside a binding site
    from variants acting somewhere else, which for effector engineering is the
    difference between a steric explanation and a conformational one.
    """
    vs = list(variants) if variants is not None else load_variants()
    cx = complexes or load_complexes()
    fcgr = [p for p, c in cx.items() if "Fcgamma" in c.partner]
    fcrn = [p for p, c in cx.items() if "FcRn" in c.partner]
    out = []
    for v in vs:
        pos = v.positions
        in_fcgr = sorted({p for p in pos if any(cx[k].at_interface(p) for k in fcgr)})
        in_fcrn = sorted({p for p in pos if any(cx[k].at_interface(p) for k in fcrn)})
        out.append(
            {
                "id": v.id,
                "label": v.label,
                "intent": v.intent,
                "isotype": v.isotype,
                "mutations": v.mutation_string(mark_disputed=True),
                "n_positions": len(pos),
                "fcgr_interface": in_fcgr,
                "fcrn_interface": in_fcrn,
                "off_interface": sorted(set(pos) - set(in_fcgr) - set(in_fcrn)),
            }
        )
    return out


# --------------------------------------------------------------------- PyMOL
_PML_HEADER = """\
# {title}
# Generated by fcatlas (Fc Variant Atlas). Reproducible: no timestamps.
#
# Structure : {pdb} - {partner}
# Fc chains : {fc}
# Fc genotype in this entry: {genotype}
{caveat}#
# Run with:  pymol {stem}.pml
# PyMOL fetches the structure over the network on first run.

reinitialize
set fetch_type_default, cif
set assembly, 1
fetch {pdb}, async=0
hide everything
bg_color white
set cartoon_transparency, 0.15
set ray_opaque_background, 0
"""


def pymol_script(
    variant: Variant,
    pdb: str = "1E4K",
    complexes: Optional[Dict[str, Complex]] = None,
    stem: Optional[str] = None,
) -> str:
    """Build a self-contained PyMOL script for one variant on one structure.

    The script fetches the structure, shows the Fc and its partner, paints the
    measured interface, then highlights and labels the variant positions. It
    also prints, in the PyMOL log, which positions are absent from the model,
    because a residue that was not resolved is not a residue that is not there.
    """
    cx = complexes or load_complexes()
    if pdb not in cx:
        raise KeyError(f"no committed interface data for {pdb}; have {sorted(cx)}")
    c = cx[pdb]
    stem = stem or f"{variant.id}_{pdb}"
    fc_sel = "+".join(c.fc_chains)
    partner_sel = "+".join(c.partner_chains)
    caveat = f"# CAVEAT  : {c.caveat}\n" if c.caveat else ""
    lines = [
        _PML_HEADER.format(
            title=f"{variant.label} ({variant.mutation_string()}) on {pdb}",
            pdb=pdb,
            partner=c.partner,
            fc=fc_sel,
            genotype=c.fc_genotype,
            caveat=caveat,
            stem=stem,
        )
    ]
    a = lines.append
    a(f'select fc, {pdb} and chain {fc_sel} and polymer')
    a(f'select partner, {pdb} and chain {partner_sel} and polymer')
    a("show cartoon, fc or partner")
    a(f"color {FC_COLOR}, fc")
    a(f"color {PARTNER_COLOR}, partner")

    iface = "+".join(str(p) for p in c.interface_positions)
    a(f"\n# measured interface, {c.cutoff} A heavy-atom cutoff")
    a(f"select interface, fc and resi {iface}")
    a(f"color {INTERFACE_COLOR}, interface")
    a("show sticks, interface and not (name C+N+O)")
    a("set stick_radius, 0.15, interface")

    pos = variant.positions
    if pos:
        sel = "+".join(str(p) for p in pos)
        a(f"\n# {variant.label}: {variant.mutation_string()}")
        a(f"select variant, fc and resi {sel}")
        a(f"color {VARIANT_COLOR}, variant")
        a("show spheres, variant and not hydrogens")
        a("set sphere_scale, 0.30, variant")
        a("show sticks, variant")
        a("set stick_radius, 0.25, variant")
        a("label variant and name CA, '%s%s' % (one_letter[resn], resi)")
        a("set label_size, 16")
        a("set label_color, black")
        a("set label_outline_color, white")

        overlap = [p for p in pos if c.at_interface(p)]
        if overlap:
            osel = "+".join(str(p) for p in overlap)
            a(f"\n# variant positions inside the measured interface: {osel}")
            a(f"select variant_at_interface, fc and resi {osel}")
            a("distance contacts, variant_at_interface, partner, 5.0, mode=0")
            a("color black, contacts")
            a("set dash_gap, 0.3")
            a("hide labels, contacts")
        else:
            a("\n# none of this variant's positions contact the partner in this")
            a("# structure, so any effect it has is not a direct steric block here")

        a("\n# report unresolved positions to the log")
        a("python")
        a(f"wanted = {pos!r}")
        a("present = set()")
        a("cmd.iterate('fc and name CA', 'present.add(int(resi))',"
          " space={'present': present, 'int': int})")
        a("missing = [p for p in wanted if p not in present]")
        a("print('fcatlas: variant positions not resolved in this model:',"
          " missing if missing else 'none')")
        a("python end")

        a(f"\norient variant or (interface within 12 of variant)")
    else:
        a(f"\n# {variant.label} carries no substitutions: {variant.notes[:70]}")
        a("orient interface")

    a("\nset antialias, 2")
    a("set ray_shadows, 0")
    a("deselect")
    return "\n".join(lines) + "\n"


def write_pymol_scripts(
    outdir: str,
    variants: Optional[Sequence[Variant]] = None,
    pdb: str = "1E4K",
) -> List[str]:
    """Write one PyMOL script per variant. Returns the paths written."""
    import os

    os.makedirs(outdir, exist_ok=True)
    vs = list(variants) if variants is not None else load_variants()
    written = []
    for v in vs:
        stem = f"{v.id}_{pdb}"
        path = os.path.join(outdir, f"{stem}.pml")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(pymol_script(v, pdb=pdb, stem=stem))
        written.append(path)
    return written
