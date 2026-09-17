# The measured interface

Every distance and area in this repository was computed from deposited
coordinates by `tools/measure_interfaces.py`. Nothing here is copied from a
paper's figure, a review table, or another database. If a number is wrong, the
script that produced it is in the repository and you can rerun it.

    python tools/measure_interfaces.py --out data/interface_detail.json

The script downloads each mmCIF from RCSB into `tools/.cache/` on first run and
reuses it afterwards. Rerunning it on unchanged inputs reproduces
`data/interface_detail.json` byte for byte.

## Method

Solvent-accessible surface area is computed with Shrake-Rupley, probe radius
1.40 A, 960 test points per atom, heavy atoms only. Waters and ions are removed
before the calculation. For each Fc chain the area is computed twice: once for
the isolated Fc chains, once for the whole complex. The difference, dSASA, is
the surface that chain gives up when the partner binds.

Contacts are heavy-atom pairs within 5.0 A. Residue numbering is the deposited
`auth_seq_id`, which is EU numbering in all four entries used here, so no
renumbering step stands between the coordinates and the table. Positions with a
dSASA below 0.5 A^2 are not reported.

Glycan chains carry their own chain identifiers in these entries. Each is
adopted by the nearest protein chain so that a sugar branching off the Fc counts
as Fc rather than as a separate binding partner, and the glycan contribution is
reported separately in `fc_glycan_contacts`.

## Two different quantities

`data/structures.json` reports, for each EU position, the minimum heavy-atom
distance to any partner atom taken over both Fc chains. It answers "is this
position anywhere near the partner in this structure".

`data/interface_detail.json` reports each position on one heavy chain: the one
that buries the most surface. Its distance and its area therefore describe the
same physical copy of the residue. It answers "how much surface does this
position give up, and to what".

The two tables can disagree, and the disagreement is informative rather than an
error. Binding to the Fc dimer is asymmetric: EU 234 buries 52.3 A^2 on one
heavy chain of 1E4K and 1.6 A^2 on the other, while EU 329 does the reverse,
120.5 A^2 against 5.2 A^2. `PositionDetail.min_distance_any_chain` carries the
cross-chain minimum so the two tables can be compared directly, and
`PositionDetail.is_asymmetric` flags positions where the two chains differ by
more than a factor of four.

## What the areas say

Proximity and burial are not the same measurement, and the atlas exists partly
to keep them apart. EU 234 sits closer to FcgammaRIIIb than EU 329 does, yet
buries a third as much surface; EU 234 also carries 17 curated variant records
against EU 329's three. Attention has followed distance, not area.

The four complexes bury 722.8 A^2 (1E4K, FcgammaRIIIb), 644.1 A^2 (1T89,
FcgammaRIIIb), 852.5 A^2 (5XJE, engineered FcgammaRIIIa) and 690.9 A^2 (4N0U,
FcRn with albumin in the asymmetric unit). Deposited glycan accounts for roughly
a tenth of each FcgammaR interface and none of the FcRn interface in 4N0U.

Two absences are worth checking against your own priors. EU 253 gives up 121.2
A^2, 96 percent of its free surface, to FcRn in 4N0U and appears in no curated
record in this atlas. EU 331 carries five records and has neither burial nor a
contact within 5 A in any of these four structures, which is what a
complement-directed position should look like in a set of structures that
contains no C1q.

## Using it from Python

    from fcatlas import load_interface, burial_by_position, engineering_against_burial

    detail = load_interface()["1E4K"]
    detail.total_buried_area          # 722.8
    detail.glycan_buried_area         # 75.0
    detail.ranked()[0].summary()      # most-buried position, one line

    p = detail.positions[329]
    p.delta_sasa, p.buried_fraction   # 120.5, 0.87
    p.is_asymmetric                   # True
    p.protein_partners[0]             # TRP87 chain C at 3.26 A

    detail.for_variant("LALA-PG")     # positions that variant changes, most buried first
    detail.variant_buried_area("LALA-PG") # 261.4

    burial_by_position()              # EU -> {pdb: area} across all complexes
    engineering_against_burial()      # one row per position: records vs area

`engineering_against_burial()` is the table behind figure 5. It joins the
curated record count for each position to the area that position buries, and it
is deliberately the plainest possible join: no ratio, no score, no ranking that
would hide which of the two numbers is doing the work.

## Scope

Four structures, resolutions 2.4 to 3.8 A. One receptor allotype per entry. The
FcgammaRIIIa complex, 5XJE, uses a receptor engineered at four positions
(N56Q, N92Q, F176V, N187Q), so it is not a wild-type FcgammaRIIIa measurement.
There is no C1q complex, no FcgammaRI, no FcgammaRIIa or IIb, and no
afucosylated Fc. An area of zero in this table means "not buried in these four
models", never "not buried".
