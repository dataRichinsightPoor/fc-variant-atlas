# Fc Variant Atlas plugin for PyMOL

Fetch a deposited Fc complex, mark the positions a named variant changes,
colour the interface by measured buried surface area, and print what the
picture rests on before it is drawn.

The plugin carries the same data as the rest of this repository: 73 curated
variant records in EU numbering, four deposited complexes with committed
contact distances, and per-position buried surface area computed from the
deposited coordinates. It needs no network access except the structure fetch
itself, and no Python packages beyond PyMOL.

## Install

Download `dist/fcatlas_pymol-<version>.zip` from the repository, then in
PyMOL choose Plugin, Plugin Manager, Install New Plugin, Choose file, and
select the zip. Restart is not required.

To run it from a clone instead, without installing:

    cd fc-variant-atlas
    pymol -r pymol_plugin/fcatlas_pymol/__init__.py

Building the zip yourself:

    python pymol_plugin/build_plugin_zip.py

The build copies `data/atlas.json`, `data/interface_detail.json` and
`data/structures.json` into the package so the installed plugin is
self-contained. If you keep the atlas installed as a Python package, or you
point `FCATLAS_DATA` at a data directory, those take precedence over the
bundled copies, so a working tree can be edited and reloaded without
rebuilding.

## Commands

    fcatlas                        the command list, data directory and version
    fcatlas_list [intent [, isotype]]
    fcatlas_show ID [, pdb [, fetch [, keep [, surface]]]]
    fcatlas_interface ID [, pdb]
    fcatlas_partners ID [, pdb [, cutoff]]
    fcatlas_burial [pdb [, object]]
    fcatlas_compare ID1, ID2 [, pdb [, offset]]
    fcatlas_provenance [pdb]
    fcatlas_hotspots [limit]

A menu entry, Fc Variant Atlas, opens a searchable variant browser on any
PyMOL built with Qt. Without Qt the commands still work; the menu is the only
thing that is missing.

## A first session

    fcatlas_show LALA-PG
    fcatlas_burial
    fcatlas_partners LALA-PG

The first command fetches 1E4K, draws the Fc and the receptor, puts the
measured interface in pale cyan and the three changed positions in red, and
prints the curated record next to the deposition's own provenance. The second
writes measured buried area into the B-factor column and colours by it, so
`select deep, b > 50` becomes a selection on buried surface rather than on
distance. The third selects and names the receptor residues within five
angstroms of the changed positions.

For a half-life variant, the default structure changes:

    fcatlas_show YTE
    fcatlas_interface YTE

The report states that the Fc chain of 4N0U carries YTE itself, that the entry
is 3.8 angstroms, that serum albumin sits in the same entry 24 angstroms from
the Fc, and that one of the three YTE substitutions, T256E, has no FcRn atom
within five angstroms and gives up 15 square angstroms of surface.

## Three choices worth knowing about

Sugars are drawn as spheres, never as lines. Glycan bond connectivity in a
deposited entry is inferred by the viewer rather than deposited, and a line
drawing asserts bonds that are not in the coordinate file. The Fc glycan is
part of the interface: it buries 75 square angstroms in 1E4K and 89 in 5XJE.

Buried area goes into the B-factor column and nothing else does, so any
selection or colouring that reads B reads measured buried area.

Burial is reported for the heavy chain that buries the most surface, with the
distance from that same chain, so one row describes one copy of the residue.
Receptor binding to the Fc homodimer is asymmetric and averaging the two
chains would invent burial on the chain that never touched the receptor. EU
234 buries 52 square angstroms on one chain and 2 on the other; EU 329 does
the reverse.

## Licence

MIT, as the rest of the repository. Structures are fetched from the PDB at
run time and are not redistributed here.
