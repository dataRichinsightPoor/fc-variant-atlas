# Open-license tools

What this repository uses, what it could be swapped for, and what has to be
cited rather than vendored. Licenses were checked against each project's own
repository or documentation.

## Structure display

**3Dmol.js** — BSD-3-Clause, [3dmol/3Dmol.js](https://github.com/3dmol/3Dmol.js).
A single script tag, no build step, WebGL rendering, works from a static file
served anywhere. This is what `viewer/index.html` uses, loaded from a pinned CDN
version. Its Python binding, `py3Dmol`, embeds the same viewer in Jupyter, which
makes a notebook version of the viewer a small piece of work rather than a port.

A note learned the hard way: 3Dmol applies an element-colored line style to
every atom on load, including HETATM groups whose bond connectivity is inferred
rather than deposited. On an Fc entry with resolved glycans this produces a dense
web of spurious bonds that looks like structure. The viewer clears all styles
first, then styles protein and glycan deliberately, and draws glycans as spheres
rather than sticks because branched-sugar connectivity inference is not reliable
enough to draw bonds from.

**Mol\*** — MIT, [molstar/molstar](https://github.com/molstar/molstar). The
engine behind the RCSB and PDBe web viewers, and the better choice if you need
volumetric data, large assemblies or trajectory support.
[MolViewSpec](https://molstar.org/mol-view-spec/) is the piece most relevant
here: a declarative JSON description of a scene — selections, colors, labels,
camera — that Mol\* renders directly. The PyMOL scripts this package generates
are procedural equivalents of exactly that, so a MolViewSpec emitter is a natural
second backend and would make the atlas viewer-agnostic.

**NGL and nglview** — MIT, [nglviewer/nglview](https://github.com/nglviewer/nglview).
Mature, fast on large structures, strong Jupyter integration.
[NGLVieweR](https://nvelden.github.io/NGLVieweR/) wraps it for R and Shiny, which
matters if your analysis stack is R rather than Python.

**PyMOL.** The open-source PyMOL source is available under a permissive license
and is what the generated `.pml` scripts target. The scripts themselves are plain
text with no dependency on any particular build, so they work with an
open-source, academic or commercial installation. This package generates PyMOL
scripts; it does not require PyMOL to be installed to generate them, and its
tests check script content rather than script execution.

## Antibody numbering

This package numbers the constant region only, deliberately. For variable
domains, use tools built for the problem.

**ANARCI** — BSD-3, [oxpig/ANARCI](https://github.com/oxpig/ANARCI). HMM-based
assignment of Chothia, Kabat, Martin, AHo and IMGT numbering to variable domains,
with germline assignment and species identification. The reference
implementation.

**ANARCI_vc** — [Fraternalilab/ANARCI_vc](https://github.com/Fraternalilab/ANARCI_vc),
a fork extending coverage to constant domains, which overlaps this repository's
territory from the other direction and is worth reading before duplicating work.

**anarci-wasm** — [tomjansen.github.io/anarci-wasm](https://tomjansen.github.io/anarci-wasm/).
ANARCI compiled to WebAssembly, running entirely in the browser. This is the
route to a numbering-aware web tool with no server, and it composes cleanly with
a 3Dmol.js or Mol\* front end.

**AbNumber** — MIT, [prihoda/AbNumber](https://github.com/prihoda/AbNumber). A
friendlier Python API over ANARCI with a chain object, alignment support and
position arithmetic. If you want variable-domain numbering inside a script
rather than a pipeline, start here.

## Therapeutic-antibody metadata

**Thera-SAbDab** — [OPIG's therapeutic antibody structural database](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/),
described in [Nucleic Acids Research](https://academic.oup.com/nar/article/48/D1/D383/5573951).
Sequences, targets, formats, clinical stage and structural coverage for named
therapeutic antibodies, downloadable in bulk. The natural join partner for this
atlas: Thera-SAbDab says which molecules exist, the atlas says what was done to
their Fc.

## IMGT: cite, do not vendor

IMGT is the authoritative source for engineered IGHG variant nomenclature and
publishes the material this curation was checked against, including the
[human IGHG numbering chart](https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGHGnber.html),
the [IGHG1 C1q, FcγR and FCGRT interaction positions](https://www.imgt.org/IMGTbiotechnology/IGHG_variant/IGHG1_C1q_FcyR_FCGRT.html),
and [FcVariantsExplorer](https://www.imgt.org/fcvariantsexplorer/).

IMGT data is subject to academic-use terms. This repository stores IMGT
engineered-variant codes as identifiers and cites IMGT as a source, and does not
redistribute IMGT tables. If you are building something downstream, read their
terms rather than inheriting an assumption from this repository.

## Sequence and structure sources

**UniProt** — constant-region sequences P01857 (IGHG1), P01859 (IGHG2), P01860
(IGHG3) and P01861 (IGHG4), CC BY 4.0. Shipped in `data/*.fasta`, trimmed at the
secreted terminus, with the trimming rule stated in the numbering module rather
than applied silently.

**RCSB PDB** — deposited coordinates and entry metadata, freely redistributable
with attribution to depositors. Structural data in this repository is derived
(contact lists), and each entry is cited by DOI.

## What is deliberately not here

No structure prediction, no affinity prediction, no developability scoring. Each
of those is a research program, and a tool that mixes measured contacts with
predicted ones invites the reader to stop distinguishing them. The atlas reports
what was measured and where the measurement came from. Predictions belong in a
layer above it, clearly labeled as such.
