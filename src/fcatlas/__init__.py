"""Fc Variant Atlas: aligned, sequence-validated, structure-linked Fc variants.

Quick start::

    from fcatlas import find, interface_report, alignment_block

    v = find("LALA-PG")[0]
    print(v.mutation_string(), v.phenotype)
    print(interface_report(v).summary())
    print(alignment_block(231, 245))

Every substitution in the dataset is checked against the UniProt sequence of
its parent isotype under EU numbering, and every structural claim is a
distance computed from deposited coordinates. Two things the package refuses
to hide: substitutions whose EU index depends on an alignment or allotype
choice, and the fact that no deposited human Fc-FcRn co-crystal structure
contains a wild-type Fc.
"""

from .numbering import (
    ACCESSIONS,
    DOMAINS,
    EU_START,
    ISOTYPES,
    EUNumbering,
    Site,
    align_isotypes,
    constant_region,
    domain_of,
    numbering,
)
from .interface import (
    BURIED_AREA_THRESHOLD,
    InterfaceDetail,
    PartnerContact,
    PositionDetail,
    burial_by_position,
    engineering_against_burial,
    load_interface,
    method_note,
)
from .structure import (
    Complex,
    InterfaceReport,
    classify,
    fcrn_provenance,
    interface_report,
    load_complexes,
    pymol_script,
    write_pymol_scripts,
)
from .variants import (
    INTENT_LABELS,
    INTENTS,
    Deletion,
    Substitution,
    Variant,
    by_intent,
    find,
    load_variants,
    position_index,
    validate_all,
)
from .export import (
    alignment_block,
    alignment_fasta,
    hotspot_table,
    to_csv,
    to_json,
    variant_rows,
)

__version__ = "1.1.0"

__all__ = [
    "BURIED_AREA_THRESHOLD",
    "InterfaceDetail",
    "PartnerContact",
    "PositionDetail",
    "burial_by_position",
    "engineering_against_burial",
    "load_interface",
    "method_note",
    "ACCESSIONS",
    "DOMAINS",
    "EU_START",
    "ISOTYPES",
    "INTENTS",
    "INTENT_LABELS",
    "Complex",
    "Deletion",
    "EUNumbering",
    "InterfaceReport",
    "Site",
    "Substitution",
    "Variant",
    "align_isotypes",
    "alignment_block",
    "alignment_fasta",
    "by_intent",
    "classify",
    "constant_region",
    "domain_of",
    "fcrn_provenance",
    "find",
    "hotspot_table",
    "interface_report",
    "load_complexes",
    "load_variants",
    "numbering",
    "position_index",
    "pymol_script",
    "to_csv",
    "to_json",
    "validate_all",
    "variant_rows",
    "write_pymol_scripts",
    "__version__",
]
