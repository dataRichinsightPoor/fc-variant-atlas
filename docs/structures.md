# Structures

## What is shipped

Four deposited complexes, each carrying measured Fc-side contacts rather than
asserted ones.

| PDB | partner | Fc chains | Fc residues in contact | Fc genotype | source |
| --- | --- | --- | --- | --- | --- |
| 1E4K | FcγRIIIb ectodomain | A, B | 19 | wild type | [10.1038/35018508](https://doi.org/10.1038/35018508) |
| 1T89 | FcγRIIIa ectodomain | A | 18 | wild type | [10.1074/jbc.M100350200](https://doi.org/10.1074/jbc.M100350200) |
| 5XJE | FcγRIIIa, both partners glycosylated | A, B | 23 | wild type | [10.1038/s41598-017-13845-8](https://doi.org/10.1038/s41598-017-13845-8) |
| 4N0U | FcRn heavy chain with β2-microglobulin | E | 16 | M252Y/S254T/T256E | [10.1074/jbc.M113.537563](https://doi.org/10.1074/jbc.M113.537563) |

## How contacts are computed

An Fc residue is in contact when any of its heavy atoms lies within 5.0 Å of any
heavy atom of any partner chain. Hydrogens are ignored because most of these
entries do not have them. Waters, ions and glycans are excluded from the partner
side; deposited glycans are kept for display but never counted as contact.

The computation buckets partner atoms into a grid of 5.0 Å cells and scans the
twenty-seven neighboring cells for each Fc atom, which makes the whole survey
cheap enough to rerun rather than trust. Residue identity and numbering come from
`auth_seq_id` and `auth_asym_id` in the mmCIF, not from the label fields, because
author numbering in these entries is EU numbering and the label numbering is not.

That equivalence was checked rather than assumed: the author-numbered residues of
1E4K reproduce the EU-numbered IgG1 sequence from the numbering module across the
full CH2 and CH3 range.

Results live in `data/structures.json` and are shipped inside the package, so
`load_complexes()` requires no network access and no structure files.

## What the glycosylated entry adds

5XJE resolves five contact positions that the deglycosylated FcγRIIIa entry 1T89
does not: EU 268, 294, 295, 296 and 326. Interface extent is not a property of
the protein pair alone; it is a property of the pair as prepared. When a variant
at EU 296 looks non-interfacial in one entry and interfacial in another, the
difference is the glycan, not the biology.

The atlas therefore reports contact per structure rather than as a single
boolean, and `fcatlas interface` prints every structure in which a position is
within cutoff, with its distance.

## The FcRn provenance problem

4N0U is the FcRn complex in this set, and its Fc is not wild type. It carries
M252Y/S254T/T256E — YTE — and all three substitutions sit inside the FcRn
contact patch. Contact distances at EU 252, 254 and 256 are therefore distances
to engineered side chains.

That prompted a survey of every human FcRn-containing entry returned by an RCSB
text search, screened for an IgG Fc chain and genotyped at the FcRn contact
positions. Composition and title of each entry were confirmed against the RCSB
entry and polymer-entity records. Twenty-one entries; results in
`data/structures.json` under `fcrn_provenance_survey`, and reachable from Python
as `fcrn_provenance()`.

Nine entries contain FcRn. Four of those nine contain an IgG Fc chain:

- 4N0U — FcRn, β2m, albumin and IgG1 Fc carrying YTE
- 7Q15 — FcRn, β2m and IgG1 Fc carrying YTE plus H433K/N434F
- 6WOL — FcRn, β2m and a monomeric IgG4 Fc carrying seven substitutions
- 6WNA — FcRn, β2m and the same monomeric IgG4 Fc lineage

Every one of those four Fc chains is engineered at FcRn contact positions. None
is wild type. The remaining five FcRn entries are serum-albumin complexes (4N0F,
4K71, 6QIO, 6QIP) or an astrovirus capsid spike complex (9NAV), with no Fc at
all.

Wild-type human Fc does appear in this survey, but unbound: 4WI2, the wild-type
member of the alanine-scanning series that mapped the FcRn site, plus 6BZ4 and
the protein A complex 1L6X. The rest of that series — 4WI3 through 4WI9, carrying
I253A, S254A, H310A, N434A, H435A, Y436A and I253A/H310A — are also unbound Fc,
despite titles that read like binding studies. They report the structural
consequence of each substitution, not the complex.

The consequence is worth stating plainly. The structural account of how a
wild-type IgG1 Fc engages FcRn is assembled from engineered Fc in complex plus
wild-type Fc alone, and the engineering sits inside the interface being
described. That does not make the account wrong. It makes it an inference with a
specific, locatable dependency, and a variant analysis that treats 4N0U as a
neutral reference frame inherits that dependency without recording it.

This claim is scoped to the survey shipped here. It is a text search plus a
composition screen, not a proof about the entire PDB, and the survey is stored
with its own `_about`, `_finding` and `_counts` fields so the scope travels with
the data.

## What the generated scripts do about it

Every PyMOL script the package emits prints the Fc genotype of the entry it
fetches. When the genotype is not wild type, the script prints a caveat block
naming the substitutions and saying that distances at those positions are
distances to engineered side chains. When a variant has no position within cutoff
in the chosen structure, the script says so in words, because an empty selection
looks identical to a negative result.

```
# CAVEAT: the Fc in this entry is not wild type (M252Y/S254T/T256E).
```

## Adding a structure

`data/structures.json` holds one object per complex with the PDB code, partner
description, chain assignments, cutoff, per-residue contacts, DOI and an optional
caveat. Contacts are generated from the mmCIF rather than hand-entered. A new
entry needs its Fc genotype determined the same way — by comparing the deposited
Fc sequence to the EU-numbered reference at every contact position — because that
genotype is what the caveat depends on.
