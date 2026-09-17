# Numbering

Three systems are in active use for the IgG constant region, they disagree, and
the disagreement changes which residue you mutate.

## EU (Edelman)

EU numbering comes from the 1969 Edelman sequence of the myeloma protein Eu and
is what the therapeutic Fc literature uses. LALA means L234A/L235A. YTE means
M252Y/S254T/T256E. DLE means S239D/A330L/I332E. When a paper says "N297A" it
means EU 297, the N-glycosylation sequon asparagine.

In this package EU numbering is defined by an alignment to human IgG1 (UniProt
P01857), trimmed at the secreted terminus, and spanning EU 118 through EU 447.
That is 330 residues for IgG1, and every one of them has an EU index.

## IMGT unique numbering

IMGT numbers the constant domain by a fixed structural framework rather than by
sequence position, so the index survives insertions and deletions across
isotypes and species. Structural and immunogenetics databases use it, and IMGT
publishes engineered-variant codes such as G1v14-49 for LALA-PG. Each record in
the atlas carries its IMGT code when one exists.

## Sequential

Sequential numbering counts residues from the first amino acid of the construct
you are actually expressing. It depends on your signal peptide, your variable
domain, your linker and your tag, so it is never portable and never appears in
this package. It is, however, what your primer design tool uses, which is why
translation errors happen at the bench rather than in the paper.

## How the atlas resolves them

The four human IgG constant regions (P01857, P01859, P01860, P01861) are trimmed
and globally aligned to IgG1 with an affine-gap Needleman–Wunsch alignment,
BLOSUM62, gap open −11, gap extend −1. The EU index is then carried as a
property of the alignment: each isotype residue either maps to an IgG1 column
with an EU number, or it does not map at all.

That last case is real. IgG3 has an extended hinge, so 47 of its 377 constant
residues have no EU counterpart in an IgG1-referenced alignment. The package
exposes them as `unnumbered_indices` rather than pretending they are numbered.

Thirty-three landmark positions — the glycosylation sequon, the hinge cysteines,
the FcRn histidines, the CH3 knob-into-hole positions and others — are asserted
against this alignment in the test suite. All thirty-three pass. The alignment is
not trusted because it is an alignment; it is trusted because it reproduces
positions whose identity is known independently.

## The IgG2 lower hinge

IgG2's lower hinge is one residue shorter than IgG1's. Under a strict alignment,
IgG2 has no residue at EU 234, and the valine that sits immediately before the
alanine is at EU 235.

```
EU     228 229 230 231 232 233 234 235 236 237
IgG1    P   C   P   A   P   E   L   L   G   G
IgG2    P   C   P   A   P   P   -   V   A   G
```

The consequence is a genuine, unresolvable-by-fiat disagreement in the
literature. Papers describing silenced IgG2 backbones write the substitution as
V234A, following EU convention as applied by the authors to their own construct.
IMGT writes the same physical residue as V235A, following the alignment. Both
statements refer to the same atom.

The atlas stores this record with `disputed: true`, reports both codes, and puts
the disputed substitution in a separate `disputed_mutations` column in the CSV
export so that a downstream join cannot silently inherit one convention. The
validator reports it as numbering-dependent rather than as an error.

`fcatlas validate` prints two flagged substitutions. The second is the E356K
member of a charge-pair heterodimerization set, where the reference residue
depends on IgG1 allotype rather than on numbering.

## Two errors this caught

An early draft of the curation in this repository contained two hand-typed
reference residues that were simply wrong.

EU 429 was recorded as proline. It is histidine. This invalidated a first pass
of the FcRn structure survey, which had been screening genotypes against the
wrong reference residue at that position.

EU 407 was recorded as threonine. It is tyrosine, so the correct notation for
that CH3 substitution is Y407V.

Neither was found by rereading. Both were found by asserting the reference
residue against the sequence at load time and letting the assertion fail. Every
reference residue in this package is now derived from the numbering module, and
the test suite forbids the literal alternative.

## Practical guidance

When you receive a variant designation, ask which system it is in before you
order the DNA. If the isotype is IgG2 or IgG3 and the position is in the lower
hinge, assume nothing. Check the residue identity, not the number:

```python
from fcatlas import numbering
numbering("IgG2").residue(234)   # None
numbering("IgG2").residue(235)   # 'V'
numbering("IgG1").residue(235)   # 'L'
```

A designation that names the wild-type residue, as LALA and YTE both do, is
self-checking. A designation that names only a position is not.
