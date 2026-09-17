"""Buried surface area and named partner residues, per Fc position.

A contact distance says a residue is near the partner. It does not say how
much of that residue the partner covers, and the two questions come apart:
in the FcgammaRIIIb complex used here, EU 234 sits 3.15 angstroms from the
receptor and gives up a third of its exposed surface, while EU 329 sits
slightly further away and gives up nearly all of it.

This module reads the committed measurements in data/interface_detail.json:
solvent-accessible surface area for every Fc residue alone and in the complex,
the difference, and every partner residue within the cutoff with its name and
distance. Areas come from Shrake-Rupley over the deposited coordinates with a
1.40 angstrom probe, heavy atoms only, glycans counted as part of whichever
chain they belong to.

Receptor binding to the Fc homodimer is asymmetric, so each position is
reported for the heavy chain that buries the most surface, with every chain's
value kept under ``per_chain``. Averaging the two chains would halve the
burial of a real interface and invent burial on a chain that never touched
the partner.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib import resources
from typing import Dict, List, Optional, Sequence, Union

from .variants import Variant, find, load_variants

#: A position giving up at least this much area is called buried.
BURIED_AREA_THRESHOLD = 10.0


@dataclass(frozen=True)
class PartnerContact:
    """One partner residue within the cutoff of an Fc position."""

    chain: str
    seq: int
    residue: str
    distance: float
    is_glycan: bool = False

    def __str__(self) -> str:
        tag = " (glycan)" if self.is_glycan else ""
        return f"{self.residue}{self.seq}[{self.chain}] {self.distance:.2f} A{tag}"


@dataclass(frozen=True)
class PositionDetail:
    """What one Fc position does at one interface."""

    eu: int
    residue: str
    chain: str
    sasa_free: float
    sasa_bound: float
    delta_sasa: float
    buried_fraction: float
    min_distance: Optional[float]
    in_contact: bool
    contacts_glycan_only: bool
    partners: List[PartnerContact] = field(default_factory=list)
    per_chain: Dict[str, dict] = field(default_factory=dict)

    @property
    def is_buried(self) -> bool:
        return self.delta_sasa >= BURIED_AREA_THRESHOLD

    @property
    def is_asymmetric(self) -> bool:
        """True when the two heavy chains bury very different amounts."""
        values = [c["delta_sasa"] for c in self.per_chain.values()]
        if len(values) < 2:
            return False
        low, high = min(values), max(values)
        return high >= BURIED_AREA_THRESHOLD and low <= 0.25 * high

    @property
    def min_distance_any_chain(self) -> Optional[float]:
        """Closest approach over both heavy chains.

        ``min_distance`` belongs to the chain that buries the most surface, so
        that the distance and the area describe the same copy of the residue.
        The committed contact table in structures.json instead takes the
        minimum over both chains, and the two disagree wherever one chain is
        closer while the other is more buried. This property is the one to
        compare against that table.
        """
        candidates = [
            c["min_distance"]
            for c in self.per_chain.values()
            if c.get("min_distance") is not None
        ]
        if self.min_distance is not None:
            candidates.append(self.min_distance)
        return min(candidates) if candidates else None

    @property
    def protein_partners(self) -> List[PartnerContact]:
        return [p for p in self.partners if not p.is_glycan]

    def summary(self) -> str:
        distance = (
            f"{self.min_distance:.2f} A" if self.min_distance is not None else "no contact"
        )
        return (
            f"EU{self.eu} {self.residue} chain {self.chain}: {distance}, "
            f"buries {self.delta_sasa:.1f} A^2 "
            f"({self.buried_fraction * 100:.0f}% of its free surface)"
        )


@dataclass(frozen=True)
class InterfaceDetail:
    """Every measured Fc position at one deposited interface."""

    pdb: str
    method: Optional[str]
    resolution_angstrom: Optional[float]
    fc_chains: List[str]
    fc_glycan_chains: List[str]
    partner_chains: List[str]
    partner_glycan_chains: List[str]
    deposited_fc_mutation: Optional[str]
    deposited_partner_mutation: Optional[str]
    positions: Dict[int, PositionDetail]
    fc_glycan_contacts: List[dict] = field(default_factory=list)
    ignored_chains: List[dict] = field(default_factory=list)

    @property
    def contact_positions(self) -> List[int]:
        return sorted(eu for eu, p in self.positions.items() if p.in_contact)

    @property
    def buried_positions(self) -> List[int]:
        return sorted(eu for eu, p in self.positions.items() if p.is_buried)

    @property
    def total_buried_area(self) -> float:
        """Fc surface the partner covers, protein plus Fc glycan."""
        protein = sum(p.delta_sasa for p in self.positions.values())
        glycan = sum(g.get("delta_sasa", 0.0) for g in self.fc_glycan_contacts)
        return round(protein + glycan, 1)

    @property
    def glycan_buried_area(self) -> float:
        return round(sum(g.get("delta_sasa", 0.0) for g in self.fc_glycan_contacts), 1)

    @property
    def wild_type_partner(self) -> bool:
        return self.deposited_partner_mutation is None

    def get(self, eu: int) -> Optional[PositionDetail]:
        return self.positions.get(eu)

    def ranked(self, limit: Optional[int] = None) -> List[PositionDetail]:
        """Positions ordered by how much surface the partner buries."""
        order = sorted(self.positions.values(), key=lambda p: -p.delta_sasa)
        return order[:limit] if limit else order

    def for_variant(self, variant: Union[Variant, str]) -> List[PositionDetail]:
        """Detail for the positions this variant changes, most buried first.

        Accepts a Variant or anything `find()` resolves: an id, a label or a
        literature alias. Raises KeyError on an unknown name rather than
        returning an empty list, because an empty list is what a variant with no
        measured position legitimately returns.
        """
        variant = _resolve(variant)
        found = [
            self.positions[eu] for eu in variant.positions if eu in self.positions
        ]
        return sorted(found, key=lambda p: -p.delta_sasa)

    def variant_buried_area(self, variant: Union[Variant, str]) -> float:
        """Total surface this variant's positions give up to the partner."""
        return round(sum(p.delta_sasa for p in self.for_variant(variant)), 1)


def _resolve(variant: Union[Variant, str]) -> Variant:
    if isinstance(variant, Variant):
        return variant
    hits = find(variant)
    if not hits:
        raise KeyError(f"no variant matches {variant!r}")
    if len(hits) > 1:
        exact = [v for v in hits if v.id.lower() == variant.strip().lower()]
        if not exact:
            names = ", ".join(v.id for v in hits[:6])
            raise KeyError(f"{variant!r} matches more than one record: {names}")
        hits = exact
    return hits[0]


def _payload() -> dict:
    return json.loads(
        resources.files("fcatlas")
        .joinpath("data", "interface_detail.json")
        .read_text(encoding="utf-8")
    )


def _position(eu: int, rec: dict) -> PositionDetail:
    return PositionDetail(
        eu=eu,
        residue=rec["residue"],
        chain=rec["chain"],
        sasa_free=rec["sasa_free"],
        sasa_bound=rec["sasa_bound"],
        delta_sasa=rec["delta_sasa"],
        buried_fraction=rec["buried_fraction"],
        min_distance=rec.get("min_distance"),
        in_contact=rec.get("in_contact", rec.get("min_distance") is not None),
        contacts_glycan_only=rec.get("contacts_glycan_only", False),
        partners=[
            PartnerContact(
                chain=p["chain"],
                seq=p["seq"],
                residue=p["residue"],
                distance=p["distance"],
                is_glycan=p.get("is_glycan", False),
            )
            for p in rec.get("partners", [])
        ],
        per_chain=rec.get("per_chain", {}),
    )


def load_interface() -> Dict[str, InterfaceDetail]:
    """Load the committed per-position interface measurements."""
    data = _payload()["complexes"]
    out: Dict[str, InterfaceDetail] = {}
    for pdb, rec in data.items():
        out[pdb] = InterfaceDetail(
            pdb=pdb,
            method=rec.get("method"),
            resolution_angstrom=rec.get("resolution_angstrom"),
            fc_chains=list(rec.get("fc_chains", [])),
            fc_glycan_chains=list(rec.get("fc_glycan_chains", [])),
            partner_chains=list(rec.get("partner_chains", [])),
            partner_glycan_chains=list(rec.get("partner_glycan_chains", [])),
            deposited_fc_mutation=rec.get("deposited_fc_mutation"),
            deposited_partner_mutation=rec.get("deposited_partner_mutation"),
            positions={
                int(eu): _position(int(eu), pos)
                for eu, pos in rec.get("positions", {}).items()
            },
            fc_glycan_contacts=list(rec.get("fc_glycan_contacts", [])),
            ignored_chains=list(rec.get("ignored_chains", [])),
        )
    return out


def method_note() -> dict:
    """How the areas and contacts were computed."""
    return _payload()["_method"]


def burial_by_position(
    interfaces: Optional[Dict[str, InterfaceDetail]] = None,
) -> Dict[int, dict]:
    """Largest burial seen at each EU position, and where it was seen."""
    interfaces = interfaces or load_interface()
    out: Dict[int, dict] = {}
    for pdb, detail in interfaces.items():
        for eu, position in detail.positions.items():
            best = out.get(eu)
            if best is None or position.delta_sasa > best["delta_sasa"]:
                out[eu] = {
                    "eu": eu,
                    "residue": position.residue,
                    "delta_sasa": position.delta_sasa,
                    "buried_fraction": position.buried_fraction,
                    "min_distance": position.min_distance,
                    "pdb": pdb,
                }
    return dict(sorted(out.items()))


def engineering_against_burial(
    variants: Optional[Sequence[Variant]] = None,
    interfaces: Optional[Dict[str, InterfaceDetail]] = None,
) -> List[dict]:
    """How often each position is engineered, against how much it buries.

    The two rankings disagree, which is the point of computing both: the
    positions the field mutates most are not the positions that give up the
    most surface when a receptor binds.
    """
    variants = variants if variants is not None else load_variants()
    burial = burial_by_position(interfaces)
    counts: Dict[int, int] = {}
    for variant in variants:
        for eu in set(variant.positions):
            counts[eu] = counts.get(eu, 0) + 1

    rows = []
    for eu in sorted(set(counts) | set(burial)):
        record = burial.get(eu)
        rows.append(
            {
                "eu": eu,
                "records": counts.get(eu, 0),
                "delta_sasa": record["delta_sasa"] if record else None,
                "buried_fraction": record["buried_fraction"] if record else None,
                "min_distance": record["min_distance"] if record else None,
                "residue": record["residue"] if record else None,
                "pdb": record["pdb"] if record else None,
            }
        )
    return rows
