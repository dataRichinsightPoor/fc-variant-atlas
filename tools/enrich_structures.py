"""Rewrite data/structures.json with provenance taken from the depositions.

Every field written here is read out of the mmCIF file rather than typed in:
chain assignment from the entity-to-chain mapping, partner identity from the
sequence database accession, genotype from the deposited mutation annotation,
resolution and method from the experiment records, and the citation DOI from
the deposition's own citation block. Contacts are left untouched; they were
recomputed independently and matched to the hundredth of an angstrom.

Two corrections came out of doing this. 1T89 was recorded with one Fc chain
and the receptor on chain B; the deposition has two Fc chains, A and B, with
the receptor on C. And 1T89's receptor is FcgammaRIIIb, not IIIa, which
leaves 5XJE as the only FcgammaRIIIa entry in the set.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

TOOLS = Path(__file__).resolve().parent
REPO = TOOLS.parent
CACHE = TOOLS / ".cache"

IG_ACCESSIONS = {"P01857", "P0DOX5", "P01859", "P01860", "P01861"}
IG_GENBANK = {"9857753"}

SPEC = {
    "1E4K": {"fc": ["A", "B"], "partner": ["C"]},
    "1T89": {"fc": ["A", "B"], "partner": ["C"]},
    "5XJE": {"fc": ["A", "B"], "partner": ["C"]},
    "4N0U": {"fc": ["E"], "partner": ["A", "B"]},
}

RECEPTOR_NAMES = {
    "O75015": "FcgammaRIIIb ectodomain",
    "P08637": "FcgammaRIIIa ectodomain",
    "P55899": "FcRn heavy chain",
    "P61769": "beta-2 microglobulin",
}

PARTNER_LABEL_SUFFIX = {"5XJE": ", glycans resolved on both partners"}

PARTNER_CAVEATS = {
    "5XJE": (
        "The receptor in this entry is not wild type. Three N-glycosylation "
        "sequons are removed (N56Q, N92Q, N187Q) and F176V is the 158V "
        "allotype in mature-protein numbering. This is the only FcgammaRIIIa "
        "entry used here; 1E4K and 1T89 both carry FcgammaRIIIb, which is "
        "97 percent identical but is not the receptor that mediates ADCC by "
        "natural killer cells."
    ),
}


def cif(pdb: str) -> dict:
    return MMCIF2Dict(str(CACHE / f"{pdb}.cif"))


def listed(d: dict, key: str) -> list:
    value = d.get(key)
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def read_provenance(pdb: str) -> dict:
    d = cif(pdb)
    chains_by_entity = defaultdict(set)
    for chain_id, entity_id in zip(
        listed(d, "_atom_site.auth_asym_id"), listed(d, "_atom_site.label_entity_id")
    ):
        chains_by_entity[entity_id].add(chain_id)

    accession_by_entity, dbname_by_entity = {}, {}
    for entity_id, accession, db in zip(
        listed(d, "_struct_ref.entity_id"),
        listed(d, "_struct_ref.pdbx_db_accession"),
        listed(d, "_struct_ref.db_name"),
    ):
        accession_by_entity[entity_id] = accession
        dbname_by_entity[entity_id] = db

    entities = []
    for i, entity_id in enumerate(listed(d, "_entity.id")):
        types = listed(d, "_entity.type")
        if i < len(types) and types[i] in {"water", "non-polymer"}:
            continue
        descriptions = listed(d, "_entity.pdbx_description")
        mutations = listed(d, "_entity.pdbx_mutation")
        mutation = mutations[i] if i < len(mutations) else "?"
        entities.append(
            {
                "entity_id": entity_id,
                "type": types[i] if i < len(types) else None,
                "description": " ".join(
                    (descriptions[i] if i < len(descriptions) else "?").split()
                )[:160],
                "accession": accession_by_entity.get(entity_id),
                "database": dbname_by_entity.get(entity_id),
                "deposited_mutation": (
                    None if mutation.strip() in {"?", "."} else mutation
                ),
                "chains": sorted(chains_by_entity.get(entity_id, ())),
            }
        )

    resolution = listed(d, "_refine.ls_d_res_high")
    method = listed(d, "_exptl.method")
    dois = [x for x in listed(d, "_citation.pdbx_database_id_DOI") if x not in {"?", "."}]
    return {
        "entities": entities,
        "resolution_angstrom": float(resolution[0]) if resolution else None,
        "method": method[0] if method else None,
        "deposited_citation_doi": f"https://doi.org/{dois[0]}" if dois else None,
    }


def main() -> None:
    payload = json.loads((REPO / "data" / "structures.json").read_text())
    detail = json.loads((REPO / "data" / "interface_detail.json").read_text())

    for pdb, spec in SPEC.items():
        record = payload["complexes"][pdb]
        prov = read_provenance(pdb)
        by_chain = {}
        for entity in prov["entities"]:
            for chain in entity["chains"]:
                by_chain[chain] = entity

        fc_entities = [by_chain[c] for c in spec["fc"]]
        for entity in fc_entities:
            accession = entity["accession"] or ""
            assert (
                accession in IG_ACCESSIONS or accession in IG_GENBANK
            ), f"{pdb}: chain assigned to Fc is not an immunoglobulin: {entity}"

        partner_entities = [by_chain[c] for c in spec["partner"]]
        partner_accessions = [e["accession"] for e in partner_entities]
        partner_label = " and ".join(
            RECEPTOR_NAMES.get(a, (e["description"] or "")[:60])
            for a, e in zip(partner_accessions, partner_entities)
        ) + PARTNER_LABEL_SUFFIX.get(pdb, "")

        other = [
            {
                "chain": e["chains"][0] if e["chains"] else None,
                "description": e["description"],
                "accession": e["accession"],
                "closest_approach_to_fc_angstrom": next(
                    (
                        ig["closest_approach_to_fc_angstrom"]
                        for ig in detail["complexes"][pdb]["ignored_chains"]
                        if ig["chain"] in e["chains"]
                    ),
                    None,
                ),
            }
            for e in prov["entities"]
            if e["type"] == "polymer"
            and not set(e["chains"]) & set(spec["fc"] + spec["partner"])
        ]

        if record["fc_chains"] != spec["fc"] or record["partner_chains"] != spec["partner"]:
            print(
                f"{pdb}: correcting chains, was fc={record['fc_chains']} "
                f"partner={record['partner_chains']}, now fc={spec['fc']} "
                f"partner={spec['partner']}"
            )
        if record["partner"] != partner_label:
            print(f"{pdb}: partner was {record['partner']!r}, now {partner_label!r}")

        record["fc_chains"] = spec["fc"]
        record["partner_chains"] = spec["partner"]
        record["partner"] = partner_label
        record["fc_glycan_chains"] = detail["complexes"][pdb]["fc_glycan_chains"]
        record["partner_glycan_chains"] = detail["complexes"][pdb][
            "partner_glycan_chains"
        ]
        record["resolution_angstrom"] = prov["resolution_angstrom"]
        record["method"] = prov["method"]
        record["deposited_citation_doi"] = prov["deposited_citation_doi"]
        # Carry the database name with the accession. The Fc chain of 1T89 is
        # referenced to GenBank rather than UniProt, so an unlabelled string
        # would read as a UniProt accession that does not exist.
        record["fc_accession"] = fc_entities[0]["accession"]
        record["fc_accession_database"] = fc_entities[0].get("database")
        record["partner_accessions"] = partner_accessions
        record["partner_accession_databases"] = [
            e.get("database") for e in partner_entities
        ]
        record["deposited_fc_mutation"] = fc_entities[0]["deposited_mutation"]
        record["deposited_partner_mutation"] = next(
            (e["deposited_mutation"] for e in partner_entities if e["deposited_mutation"]),
            None,
        )
        record["partner_genotype"] = (
            record["deposited_partner_mutation"].replace(", ", "/")
            if record["deposited_partner_mutation"]
            else "wild type"
        )
        record["wild_type_partner"] = record["deposited_partner_mutation"] is None
        record["partner_caveat"] = PARTNER_CAVEATS.get(pdb, "")
        record["other_polymers_in_entry"] = other
        record["entities"] = prov["entities"]

        deposited = record["deposited_fc_mutation"]
        curated = record["fc_genotype"]
        if deposited:
            for sub in deposited.replace(" ", "").split(","):
                assert sub in curated.replace(" ", ""), (
                    f"{pdb}: deposited Fc mutation {sub} is missing from the curated "
                    f"genotype {curated!r}"
                )
            print(f"{pdb}: deposited Fc mutation {deposited!r} agrees with {curated!r}")
        elif curated.strip().lower() != "wild type":
            raise AssertionError(
                f"{pdb}: curated genotype is {curated!r} but the deposition annotates "
                "no mutation"
            )

    marker = " Chain assignment"
    payload["_method"] = payload["_method"].split(marker)[0]
    payload["_method"] = (
        payload["_method"]
        + " Chain assignment, partner identity, genotype, resolution, method and "
        "citation are read from the deposited mmCIF rather than transcribed: the "
        "entity-to-chain mapping gives the chains, the sequence database accession "
        "gives the partner, and _entity.pdbx_mutation gives the genotype, which is "
        "checked against the curated genotype string at build time."
    )
    payload["_partner_provenance_finding"] = (
        "Two of the three FcgammaR complexes used here, 1E4K and 1T89, carry "
        "FcgammaRIIIb (UniProt O75015). Only 5XJE carries FcgammaRIIIa (P08637), "
        "the receptor responsible for natural killer cell ADCC, and its receptor is "
        "engineered: N56Q, N92Q and N187Q remove three N-glycosylation sequons and "
        "F176V is the 158V allotype. The provenance problem is therefore symmetric. "
        "It is not only the Fc in these complexes that is engineered."
    )

    out = REPO / "data" / "structures.json"
    out.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {out}")

    target = REPO / "data" / "interface_detail.json"
    target.write_text(json.dumps(detail, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
