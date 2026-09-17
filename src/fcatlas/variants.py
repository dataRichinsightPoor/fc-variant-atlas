"""Curated Fc variant records, parsed and checked against real sequence.

The central guarantee of this module is that every substitution in the
dataset is validated against the actual UniProt sequence of its parent
isotype under EU numbering. A record whose stated wild-type residue does not
match the sequence is a data error, and `validate_all` finds it rather than
letting it propagate into a figure or a structure selection.

Two record classes are deliberately allowed to carry no substitutions:
glycan variants such as afucosylation, and designs containing residue
deletions. They are flagged, not silently dropped, because their existence
is the main limitation of mutation-set thinking about Fc engineering.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from importlib import resources
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .numbering import numbering

AA = set("ACDEFGHIKLMNPQRSTVWY")

_SUB_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
_DEL_RE = re.compile(r"^([A-Z])(\d+)del$")

INTENTS = (
    "silence",
    "silence+stabilize",
    "silence+halflife",
    "enhance_adcc",
    "enhance_cdc",
    "reduce_cdc",
    "enhance_fcgr2b",
    "halflife",
    "reduce_halflife",
    "heterodimer",
    "conjugation",
    "stabilize",
)

#: Human-readable labels for the engineering intents.
INTENT_LABELS = {
    "silence": "effector silencing",
    "silence+stabilize": "silencing with hinge stabilization",
    "silence+halflife": "silencing with half-life extension",
    "enhance_adcc": "ADCC / ADCP enhancement",
    "enhance_cdc": "complement enhancement",
    "reduce_cdc": "complement reduction",
    "enhance_fcgr2b": "inhibitory FcgammaRIIb enhancement",
    "halflife": "half-life extension",
    "reduce_halflife": "half-life reduction",
    "heterodimer": "heavy-chain heterodimerization",
    "conjugation": "site-specific conjugation",
    "stabilize": "stabilization",
}


@dataclass(frozen=True)
class Substitution:
    """One substitution in EU numbering."""

    wt: str
    eu: int
    mut: str

    @classmethod
    def parse(cls, text: str) -> "Substitution":
        m = _SUB_RE.match(text)
        if not m:
            raise ValueError(f"not a substitution: {text!r}")
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        if wt not in AA or mut not in AA:
            raise ValueError(f"non-standard residue in {text!r}")
        return cls(wt, pos, mut)

    def __str__(self) -> str:
        return f"{self.wt}{self.eu}{self.mut}"


@dataclass(frozen=True)
class Deletion:
    """A residue deletion in EU numbering."""

    wt: str
    eu: int

    @classmethod
    def parse(cls, text: str) -> "Deletion":
        m = _DEL_RE.match(text)
        if not m:
            raise ValueError(f"not a deletion: {text!r}")
        return cls(m.group(1), int(m.group(2)))

    def __str__(self) -> str:
        return f"{self.wt}{self.eu}del"


@dataclass
class Variant:
    """One curated Fc variant record."""

    id: str
    isotype: str
    substitutions: List[Substitution] = field(default_factory=list)
    deletions: List[Deletion] = field(default_factory=list)
    #: Substitutions that are real but whose wild-type residue or EU index is
    #: not reproducible against this package's reference frame. Two causes
    #: occur in practice: a numbering-convention collision in an isotype whose
    #: hinge cannot be aligned to IgG1 without a gap, and an allotype
    #: difference between the published construct and the UniProt reference.
    #: They are carried, displayed and counted, but exempt from strict
    #: sequence validation, and `notes` states which cause applies.
    disputed: List[Substitution] = field(default_factory=list)
    aliases: List[str] = field(default_factory=list)
    imgt: Optional[str] = None
    intent: str = ""
    phenotype: str = ""
    evidence: str = ""
    therapeutics: List[str] = field(default_factory=list)
    source: str = ""
    notes: str = ""

    # ------------------------------------------------------------- properties
    @property
    def positions(self) -> List[int]:
        """All EU positions touched, including disputed ones."""
        return sorted(
            {s.eu for s in self.substitutions}
            | {d.eu for d in self.deletions}
            | {s.eu for s in self.disputed}
        )

    @property
    def validated_positions(self) -> List[int]:
        """EU positions confirmed against the reference sequence."""
        return sorted({s.eu for s in self.substitutions})

    @property
    def is_disputed(self) -> bool:
        return bool(self.disputed)

    @property
    def is_glycan_variant(self) -> bool:
        """True for records that change glycan, not sequence."""
        return not self.substitutions and not self.deletions and not self.disputed

    @property
    def has_deletion(self) -> bool:
        return bool(self.deletions)

    @property
    def is_substitutional(self) -> bool:
        """True when the variant is fully expressible as point substitutions."""
        return bool(self.substitutions) and not self.deletions

    @property
    def label(self) -> str:
        return self.aliases[0] if self.aliases else self.id

    def mutation_string(self, sep: str = "/", mark_disputed: bool = False) -> str:
        parts = [str(s) for s in self.substitutions] + [str(d) for d in self.deletions]
        parts += [f"{s}*" if mark_disputed else str(s) for s in self.disputed]
        return sep.join(parts)

    # ------------------------------------------------------------- validation
    def validate(self) -> List[str]:
        """Return a list of problems; empty means the record is consistent."""
        problems: List[str] = []
        num = numbering(self.isotype)
        for s in self.substitutions:
            actual = num.residue(s.eu)
            if actual is None:
                problems.append(
                    f"{self.id}: EU {s.eu} has no {self.isotype} counterpart "
                    f"(record says {s.wt}); the position falls in an alignment gap"
                )
            elif actual != s.wt:
                problems.append(
                    f"{self.id}: EU {s.eu} is {actual} in {self.isotype}, "
                    f"record says {s.wt}"
                )
            if s.wt == s.mut:
                problems.append(f"{self.id}: {s} is not a substitution")
        for d in self.deletions:
            actual = num.residue(d.eu)
            if actual is not None and actual != d.wt:
                problems.append(
                    f"{self.id}: EU {d.eu} is {actual} in {self.isotype}, "
                    f"record says {d.wt}"
                )
        if self.intent not in INTENTS:
            problems.append(f"{self.id}: unknown intent {self.intent!r}")
        if not self.source:
            problems.append(f"{self.id}: missing source")
        if not self.phenotype:
            problems.append(f"{self.id}: missing phenotype")
        if self.disputed and "DISPUTE" not in self.notes.upper() and (
            "ALLOTYPE" not in self.notes.upper()
        ):
            problems.append(
                f"{self.id}: carries disputed substitutions but notes do not "
                "explain the numbering or allotype cause"
            )
        return problems

    def mutated_sequence(self) -> str:
        """Apply the substitutions to the parent sequence.

        Deletions are not applied; a variant carrying one raises, because
        silently returning a sequence of the wrong length would be worse than
        refusing.
        """
        if self.deletions:
            raise ValueError(
                f"{self.id} contains a deletion and cannot be rendered as a "
                "substituted sequence"
            )
        num = numbering(self.isotype)
        seq = list(num.sequence)
        for s in self.substitutions:
            idx = num.index(s.eu)
            if idx is None:
                raise ValueError(f"{self.id}: EU {s.eu} absent in {self.isotype}")
            seq[idx] = s.mut
        return "".join(seq)


# ------------------------------------------------------------------ minimal YAML
def _parse_yaml(text: str) -> List[Dict[str, object]]:
    """Parse the restricted YAML subset used by data/variants.yaml.

    Supports: a top-level `variants:` key, a sequence of mappings, scalar
    values, inline flow lists, and folded block scalars introduced by `>`.
    Written in-package so the dataset loads with no third-party dependency.
    """
    lines = text.splitlines()
    records: List[Dict[str, object]] = []
    cur: Optional[Dict[str, object]] = None
    i = 0
    in_variants = False
    while i < len(lines):
        raw = lines[i]
        line = raw.split("#", 1)[0].rstrip() if not raw.strip().startswith("#") else ""
        if not line.strip():
            i += 1
            continue
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())
        if stripped == "variants:":
            in_variants = True
            i += 1
            continue
        if not in_variants:
            i += 1
            continue
        if stripped.startswith("- "):
            cur = {}
            records.append(cur)
            stripped = stripped[2:]
            indent += 2
        if cur is None:
            i += 1
            continue
        if ":" not in stripped:
            i += 1
            continue
        key, _, val = stripped.partition(":")
        key = key.strip()
        val = val.strip()
        if val == ">":
            block: List[str] = []
            i += 1
            while i < len(lines):
                nxt = lines[i]
                if not nxt.strip():
                    block.append("")
                    i += 1
                    continue
                nindent = len(nxt) - len(nxt.lstrip())
                if nindent <= indent:
                    break
                block.append(nxt.strip())
                i += 1
            cur[key] = " ".join(p for p in block if p).strip()
            continue
        if val.startswith("[") and val.endswith("]"):
            inner = val[1:-1].strip()
            cur[key] = [p.strip() for p in inner.split(",") if p.strip()] if inner else []
        elif val.startswith('"') and val.endswith('"') and len(val) >= 2:
            cur[key] = val[1:-1]
        else:
            cur[key] = val
        i += 1
    return records


def _mutation_tokens(raw: object) -> Tuple[List[Substitution], List[Deletion]]:
    subs: List[Substitution] = []
    dels: List[Deletion] = []
    if not raw:
        return subs, dels
    tokens = raw if isinstance(raw, list) else [str(raw)]
    for tok in tokens:
        tok = str(tok).strip()
        if not tok:
            continue
        if tok.endswith("del"):
            dels.append(Deletion.parse(tok))
        else:
            subs.append(Substitution.parse(tok))
    return subs, dels


def load_variants(path: Optional[str] = None) -> List[Variant]:
    """Load the curated variant dataset."""
    if path is None:
        text = (
            resources.files("fcatlas")
            .joinpath("data", "variants.yaml")
            .read_text(encoding="utf-8")
        )
    else:
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
    out: List[Variant] = []
    for rec in _parse_yaml(text):
        subs, dels = _mutation_tokens(rec.get("mutations"))
        disputed, _ = _mutation_tokens(rec.get("disputed_mutations"))
        imgt = rec.get("imgt") or None
        out.append(
            Variant(
                id=str(rec["id"]),
                isotype=str(rec["isotype"]),
                substitutions=subs,
                deletions=dels,
                disputed=disputed,
                aliases=list(rec.get("aliases") or []),
                imgt=str(imgt) if imgt else None,
                intent=str(rec.get("intent", "")),
                phenotype=str(rec.get("phenotype", "")),
                evidence=str(rec.get("evidence", "")),
                therapeutics=list(rec.get("therapeutics") or []),
                source=str(rec.get("source", "")),
                notes=str(rec.get("notes", "")),
            )
        )
    return out


def validate_all(variants: Optional[Sequence[Variant]] = None) -> List[str]:
    """Validate every record against sequence. Empty result means clean."""
    vs = list(variants) if variants is not None else load_variants()
    problems: List[str] = []
    seen: Dict[str, str] = {}
    for v in vs:
        if v.id in seen:
            problems.append(f"duplicate id {v.id}")
        seen[v.id] = v.isotype
        problems.extend(v.validate())
    return problems


def by_intent(variants: Optional[Iterable[Variant]] = None) -> Dict[str, List[Variant]]:
    vs = list(variants) if variants is not None else load_variants()
    out: Dict[str, List[Variant]] = {}
    for v in vs:
        out.setdefault(v.intent, []).append(v)
    return out


def find(query: str, variants: Optional[Iterable[Variant]] = None) -> List[Variant]:
    """Look a variant up by id, alias, IMGT code, or mutation string."""
    vs = list(variants) if variants is not None else load_variants()
    q = query.strip().lower()
    hits = []
    for v in vs:
        keys = {v.id.lower(), (v.imgt or "").lower(), v.mutation_string().lower()}
        keys |= {a.lower() for a in v.aliases}
        keys.discard("")
        if q in keys:
            hits.append(v)
    if hits:
        return hits
    return [
        v
        for v in vs
        if q in v.id.lower() or any(q in a.lower() for a in v.aliases)
    ]


def position_index(
    variants: Optional[Iterable[Variant]] = None,
) -> Dict[int, List[Variant]]:
    """EU position -> variants touching it. The basis of the hotspot view."""
    vs = list(variants) if variants is not None else load_variants()
    out: Dict[int, List[Variant]] = {}
    for v in vs:
        for p in v.positions:
            out.setdefault(p, []).append(v)
    return dict(sorted(out.items()))
