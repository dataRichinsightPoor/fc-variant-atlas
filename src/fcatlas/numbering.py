"""EU numbering for human IGHG constant regions.

The EU index (Edelman 1969) is the convention used by essentially all Fc
engineering literature, but it is defined on IgG1. Applying it to IgG2, IgG3
and IgG4 requires an alignment, because the hinges differ in length. This
module builds that alignment explicitly and exposes where it is ambiguous
rather than hiding the ambiguity behind a fixed offset.

The IgG1 secreted constant region begins at EU 118 (A of ASTKGPSVFPLAP) and
ends at EU 447 (K of SLSLSPGK). That anchor is validated in the test suite
against known landmark residues and against author numbering in deposited
Fc structures, where auth_seq_id is EU numbering.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources
from typing import Dict, Iterator, List, Optional, Tuple

ISOTYPES = ("IgG1", "IgG2", "IgG3", "IgG4")

#: UniProt accessions for the human heavy-chain constant regions.
ACCESSIONS = {
    "IgG1": "P01857",
    "IgG2": "P01859",
    "IgG3": "P01860",
    "IgG4": "P01861",
}

#: First residue of the secreted constant region in EU numbering.
EU_START = 118

#: The UniProt entries continue into the membrane-anchor form. The secreted
#: allotype ends just before this motif; we trim there and restore the
#: C-terminal GK that the secreted protein carries.
_MEMBRANE_MOTIF = "ELQLEESCAEA"

# Domain boundaries in EU numbering. CH1 and the hinge are included because
# several bispecific and conjugation variants are defined outside the Fc.
DOMAINS: Tuple[Tuple[str, int, int], ...] = (
    ("CH1", 118, 215),
    ("hinge", 216, 230),
    ("CH2", 231, 340),
    ("CH3", 341, 447),
)


def domain_of(eu: int) -> str:
    """Return the domain containing an EU position."""
    for name, lo, hi in DOMAINS:
        if lo <= eu <= hi:
            return name
    return "outside"


@functools.lru_cache(maxsize=None)
def _raw_sequence(isotype: str) -> str:
    acc = ACCESSIONS[isotype]
    text = (
        resources.files("fcatlas")
        .joinpath("data", f"{acc}.fasta")
        .read_text(encoding="utf-8")
    )
    return "".join(
        line.strip() for line in text.splitlines() if not line.startswith(">")
    )


@functools.lru_cache(maxsize=None)
def constant_region(isotype: str) -> str:
    """Secreted heavy-chain constant region for an isotype."""
    if isotype not in ACCESSIONS:
        raise KeyError(f"unknown isotype {isotype!r}; expected one of {ISOTYPES}")
    raw = _raw_sequence(isotype)
    cut = raw.find(_MEMBRANE_MOTIF)
    if cut < 0:
        raise ValueError(f"membrane motif not found in {isotype}")
    return raw[:cut] + "GK"


@dataclass(frozen=True)
class Site:
    """One aligned position."""

    eu: int
    residue: str
    index: int  # zero-based index into the isotype's own sequence
    domain: str

    def __str__(self) -> str:  # pragma: no cover - display helper
        return f"{self.residue}{self.eu}"


class EUNumbering:
    """EU-numbered view of one isotype's constant region.

    For IgG1 the mapping is the defining anchor. For the other isotypes it is
    derived by global alignment to IgG1, so some EU positions have no
    counterpart (the IgG1 residue is deleted in that isotype) and some
    residues have no EU number (an insertion relative to IgG1).
    """

    def __init__(self, isotype: str) -> None:
        self.isotype = isotype
        self.sequence = constant_region(isotype)
        self._by_eu: Dict[int, Site] = {}
        self._unnumbered: List[int] = []
        if isotype == "IgG1":
            for i, aa in enumerate(self.sequence):
                eu = EU_START + i
                self._by_eu[eu] = Site(eu, aa, i, domain_of(eu))
        else:
            self._build_by_alignment()

    def _build_by_alignment(self) -> None:
        ref = EUNumbering("IgG1") if self.isotype != "IgG1" else None
        assert ref is not None
        aligned = _align(ref.sequence, self.sequence)
        i_ref = i_self = 0
        for a, b in aligned:
            if a != "-" and b != "-":
                eu = EU_START + i_ref
                self._by_eu[eu] = Site(eu, b, i_self, domain_of(eu))
            elif a == "-" and b != "-":
                self._unnumbered.append(i_self)
            if a != "-":
                i_ref += 1
            if b != "-":
                i_self += 1

    # ------------------------------------------------------------------ access
    def residue(self, eu: int) -> Optional[str]:
        """Residue at an EU position, or None if absent in this isotype."""
        site = self._by_eu.get(eu)
        return site.residue if site else None

    def site(self, eu: int) -> Optional[Site]:
        return self._by_eu.get(eu)

    def index(self, eu: int) -> Optional[int]:
        site = self._by_eu.get(eu)
        return site.index if site else None

    def positions(self) -> List[int]:
        return sorted(self._by_eu)

    @property
    def unnumbered_indices(self) -> List[int]:
        """Sequence indices that have no EU counterpart (insertions vs IgG1)."""
        return list(self._unnumbered)

    def __iter__(self) -> Iterator[Site]:
        for eu in self.positions():
            yield self._by_eu[eu]

    def __len__(self) -> int:
        return len(self._by_eu)

    def __repr__(self) -> str:  # pragma: no cover - display helper
        return (
            f"<EUNumbering {self.isotype} "
            f"{len(self._by_eu)} numbered / {len(self.sequence)} residues>"
        )


@functools.lru_cache(maxsize=None)
def numbering(isotype: str) -> EUNumbering:
    """Cached EU numbering for an isotype."""
    return EUNumbering(isotype)


# --------------------------------------------------------------------- aligner
#: BLOSUM62, embedded so the package has no hard dependency on Biopython.
_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62_ROWS = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]


@functools.lru_cache(maxsize=1)
def _blosum62() -> Dict[Tuple[str, str], int]:
    table: Dict[Tuple[str, str], int] = {}
    for i, row in enumerate(_B62_ROWS):
        for j, v in enumerate(row):
            table[(_B62_ORDER[i], _B62_ORDER[j])] = v
    return table


def _score(a: str, b: str) -> int:
    return _blosum62().get((a, b), -4)


def _align(
    ref: str, query: str, gap_open: int = -11, gap_extend: int = -1
) -> List[Tuple[str, str]]:
    """Global affine-gap alignment (Gotoh), returned as aligned column pairs.

    Implemented in-package so cross-isotype numbering is reproducible without
    an external alignment dependency or a BLAST-style heuristic.
    """
    n, m = len(ref), len(query)
    neg = float("-inf")
    # M: match state, X: gap in query, Y: gap in ref
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    X = [[neg] * (m + 1) for _ in range(n + 1)]
    Y = [[neg] * (m + 1) for _ in range(n + 1)]
    PM = [[0] * (m + 1) for _ in range(n + 1)]
    PX = [[0] * (m + 1) for _ in range(n + 1)]
    PY = [[0] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
        PX[i][0] = 1
    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend
        PY[0][j] = 2
    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(1, m + 1):
            s = _score(ri, query[j - 1])
            best, ptr = M[i - 1][j - 1], 0
            if X[i - 1][j - 1] > best:
                best, ptr = X[i - 1][j - 1], 1
            if Y[i - 1][j - 1] > best:
                best, ptr = Y[i - 1][j - 1], 2
            M[i][j], PM[i][j] = best + s, ptr

            open_m = M[i - 1][j] + gap_open
            ext_x = X[i - 1][j] + gap_extend
            if open_m >= ext_x:
                X[i][j], PX[i][j] = open_m, 0
            else:
                X[i][j], PX[i][j] = ext_x, 1

            open_m2 = M[i][j - 1] + gap_open
            ext_y = Y[i][j - 1] + gap_extend
            if open_m2 >= ext_y:
                Y[i][j], PY[i][j] = open_m2, 0
            else:
                Y[i][j], PY[i][j] = ext_y, 2

    state = max(((M[n][m], 0), (X[n][m], 1), (Y[n][m], 2)))[1]
    i, j = n, m
    out: List[Tuple[str, str]] = []
    while i > 0 or j > 0:
        if state == 0:
            out.append((ref[i - 1], query[j - 1]))
            state = PM[i][j]
            i, j = i - 1, j - 1
        elif state == 1:
            out.append((ref[i - 1], "-"))
            state = PX[i][j]
            i -= 1
        else:
            out.append(("-", query[j - 1]))
            state = PY[i][j]
            j -= 1
        if i == 0 and j > 0:
            state = 2
        elif j == 0 and i > 0:
            state = 1
    out.reverse()
    return out


def align_isotypes(
    isotypes: Tuple[str, ...] = ISOTYPES,
) -> List[Tuple[Optional[int], Dict[str, str]]]:
    """Column-wise alignment of the isotypes in EU coordinates.

    Returns a list of (eu_position, {isotype: residue_or_gap}). Positions
    present in IgG1 carry an EU number; insertions relative to IgG1 carry
    None, which is the honest representation of an unnumberable residue.
    """
    ref = numbering("IgG1")
    rows: List[Tuple[Optional[int], Dict[str, str]]] = []
    for eu in ref.positions():
        col = {}
        for iso in isotypes:
            col[iso] = numbering(iso).residue(eu) or "-"
        rows.append((eu, col))
    return rows
