"""Recompute the per-position interface measurements from deposited coordinates.

    python tools/measure_interfaces.py

Downloads the mmCIF entries it needs into tools/.cache (once), measures every
Fc residue, and writes data/interface_detail.json. Run
tools/enrich_structures.py afterwards to fold the deposition provenance into
data/structures.json. Requires biopython.

For each deposited complex, and for every Fc residue, compute:
  * solvent-accessible surface area in the isolated Fc, with its own glycans
  * the same area in the complex
  * the difference, which is the surface the partner buries
  * every partner residue within the cutoff, named, with its distance

and record the provenance the deposition itself states: resolution, method,
the polymer entities present, and the mutation annotation on each entity,
which is where the Fc genotype and the partner genotype come from.

Shrake-Rupley, probe radius 1.40 A, 960 test points, heavy atoms only,
waters and ions removed. Glycan-only chains are adopted by whichever protein
chain they sit closest to, so an Fc glycan counts as Fc surface. Polymer
chains that belong to neither side of the measured pair are excluded and
reported with their closest approach to the Fc, so their presence in the
crystal is visible rather than silently averaged in.

Each Fc position is reported for the chain that buries the most surface,
because receptor binding to the Fc homodimer is asymmetric and averaging the
two heavy chains would halve the burial of a real interface.
"""

from __future__ import annotations

import copy
import json
import math
import warnings
from collections import defaultdict

from pathlib import Path
from urllib.request import urlopen

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings("ignore")

TOOLS = Path(__file__).resolve().parent
REPO = TOOLS.parent
CACHE = TOOLS / ".cache"
RCSB = "https://files.rcsb.org/download/{}.cif"


def cif_path(pdb: str) -> str:
    """Path to a local mmCIF, downloading it from the PDB on first use."""
    CACHE.mkdir(exist_ok=True)
    target = CACHE / f"{pdb}.cif"
    if not target.is_file():
        url = RCSB.format(pdb.upper())
        print(f"fetching {url}")
        with urlopen(url, timeout=120) as response:
            target.write_bytes(response.read())
    return str(target)

CUTOFF = 5.0
PROBE = 1.40
N_POINTS = 960
MIN_DELTA = 0.5

SOLVENT_AND_IONS = {
    "HOH", "DOD", "CL", "NA", "K", "MG", "CA", "ZN", "SO4", "PO4", "GOL",
    "EDO", "ACT", "MN", "NI", "CD", "BR", "IOD", "FMT", "PEG", "TRS", "MES",
}
GLYCANS = {
    "NAG", "NDG", "BMA", "MAN", "FUC", "FUL", "GAL", "GLA", "SIA", "BGC",
    "XYS", "A2G", "GLC", "NGA",
}

JOBS = {
    "1E4K": {"fc": ["A", "B"], "partner": ["C"], "partner_name": "FcgammaRIIIb"},
    "1T89": {"fc": ["A", "B"], "partner": ["C"], "partner_name": "FcgammaRIIIa"},
    "5XJE": {"fc": ["A", "B"], "partner": ["C"], "partner_name": "FcgammaRIIIa"},
    "4N0U": {"fc": ["E"], "partner": ["A", "B"], "partner_name": "FcRn"},
}


def is_wanted(residue) -> bool:
    name = residue.get_resname().strip()
    if name in SOLVENT_AND_IONS:
        return False
    if residue.id[0] == " ":
        return True
    return name in GLYCANS


def heavy(residue):
    return [a for a in residue if a.element != "H"]


def load_model(path):
    model = MMCIFParser(QUIET=True).get_structure("x", path)[0]
    for chain in list(model):
        for residue in list(chain):
            if not is_wanted(residue):
                chain.detach_child(residue.id)
                continue
            for atom in list(residue):
                if atom.element == "H":
                    residue.detach_child(atom.id)
        if len(chain) == 0:
            model.detach_child(chain.id)
    return model


def metadata(path):
    """Entities, chains, mutation annotations and experiment details."""
    d = MMCIF2Dict(path)

    def listed(key):
        value = d.get(key)
        if value is None:
            return []
        return value if isinstance(value, list) else [value]

    chains_by_entity = defaultdict(set)
    for chain_id, entity_id in zip(
        listed("_atom_site.auth_asym_id"), listed("_atom_site.label_entity_id")
    ):
        chains_by_entity[entity_id].add(chain_id)

    entities = []
    ids = listed("_entity.id")
    descriptions = listed("_entity.pdbx_description")
    mutations = listed("_entity.pdbx_mutation")
    types = listed("_entity.type")
    for i, entity_id in enumerate(ids):
        description = descriptions[i] if i < len(descriptions) else "?"
        mutation = mutations[i] if i < len(mutations) else "?"
        entity_type = types[i] if i < len(types) else "?"
        if entity_type in {"water", "non-polymer"}:
            continue
        entities.append(
            {
                "entity_id": entity_id,
                "type": entity_type,
                "description": " ".join(description.split())[:200],
                "deposited_mutation": None if mutation.strip() in {"?", "."} else mutation,
                "chains": sorted(chains_by_entity.get(entity_id, ())),
            }
        )

    resolution = listed("_refine.ls_d_res_high") or listed(
        "_em_3d_reconstruction.resolution"
    )
    method = listed("_exptl.method")
    return {
        "method": method[0] if method else None,
        "resolution_angstrom": float(resolution[0]) if resolution else None,
        "entities": entities,
    }


def classify_chains(model, fc_protein, partner_protein):
    """Adopt glycan-only chains; leave other polymers out and report them."""
    protein_atoms = {
        chain.id: [a for r in chain if r.id[0] == " " for a in heavy(r)]
        for chain in model
        if chain.id in set(fc_protein) | set(partner_protein)
    }
    fc_glycan, partner_glycan, ignored = [], [], []
    for chain in model:
        if chain.id in fc_protein or chain.id in partner_protein:
            continue
        residues = list(chain)
        atoms = [a for r in residues for a in heavy(r)]
        if not atoms:
            continue
        glycan_only = all(r.get_resname().strip() in GLYCANS for r in residues)
        nearest_chain, nearest = None, 1e9
        for cid, patoms in protein_atoms.items():
            for a in atoms:
                for b in patoms:
                    d = a - b
                    if d < nearest:
                        nearest, nearest_chain = d, cid
        if glycan_only:
            (fc_glycan if nearest_chain in fc_protein else partner_glycan).append(
                chain.id
            )
        else:
            fc_atoms = [a for cid in fc_protein for a in protein_atoms[cid]]
            closest_to_fc = min(a - b for a in atoms for b in fc_atoms)
            ignored.append(
                {
                    "chain": chain.id,
                    "residues": len(residues),
                    "closest_approach_to_fc_angstrom": round(float(closest_to_fc), 2),
                    "reason": "polymer outside the measured pair",
                }
            )
    return sorted(fc_glycan), sorted(partner_glycan), ignored


def sasa_by_residue(model, keep_chains):
    subset = copy.deepcopy(model)
    for chain in list(subset):
        if chain.id not in keep_chains:
            subset.detach_child(chain.id)
    ShrakeRupley(probe_radius=PROBE, n_points=N_POINTS).compute(subset, level="R")
    return {
        (chain.id, residue.id[1], residue.get_resname().strip()): residue.sasa
        for chain in subset
        for residue in chain
    }


def partner_contacts(model, fc_chains, partner_chains):
    grid = defaultdict(list)
    for chain in model:
        if chain.id not in partner_chains:
            continue
        for residue in chain:
            for atom in heavy(residue):
                x, y, z = atom.coord
                grid[
                    (int(x // CUTOFF), int(y // CUTOFF), int(z // CUTOFF))
                ].append((chain.id, residue, atom))

    found = defaultdict(dict)
    for chain in model:
        if chain.id not in fc_chains:
            continue
        for residue in chain:
            key = (chain.id, residue.id[1], residue.get_resname().strip())
            for atom in heavy(residue):
                x, y, z = atom.coord
                gx, gy, gz = int(x // CUTOFF), int(y // CUTOFF), int(z // CUTOFF)
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dz in (-1, 0, 1):
                            for pchain, presidue, patom in grid.get(
                                (gx + dx, gy + dy, gz + dz), ()
                            ):
                                d = math.dist(atom.coord, patom.coord)
                                if d > CUTOFF:
                                    continue
                                pkey = (
                                    pchain,
                                    presidue.id[1],
                                    presidue.get_resname().strip(),
                                )
                                if pkey not in found[key] or d < found[key][pkey]:
                                    found[key][pkey] = d
    return found


def build_positions(free, bound, contacts, fc_protein, fc_glycan):
    """One record per EU position, from the chain that buries the most."""
    per_position = defaultdict(dict)
    glycan_records = []

    for key in sorted(free, key=lambda k: (k[1], k[0])):
        chain_id, seq_id, resname = key
        free_area = free[key]
        bound_area = bound.get(key, free_area)
        delta = free_area - bound_area
        partners = [
            {
                "chain": pchain,
                "seq": pseq,
                "residue": pname,
                "distance": round(float(d), 2),
                "is_glycan": pname in GLYCANS,
            }
            for (pchain, pseq, pname), d in sorted(
                contacts.get(key, {}).items(), key=lambda kv: kv[1]
            )
        ]
        if delta < MIN_DELTA and not partners:
            continue
        record = {
            "residue": resname,
            "chain": chain_id,
            "sasa_free": round(float(free_area), 1),
            "sasa_bound": round(float(bound_area), 1),
            "delta_sasa": round(float(delta), 1),
            "buried_fraction": (
                round(float(delta) / float(free_area), 3) if free_area > 0 else 0.0
            ),
            "min_distance": partners[0]["distance"] if partners else None,
            "partners": partners,
        }
        if chain_id in fc_glycan:
            record["glycan_chain"] = chain_id
            record["seq"] = seq_id
            glycan_records.append(record)
        else:
            per_position[seq_id][chain_id] = record

    positions = {}
    for seq_id, by_chain in per_position.items():
        representative = max(by_chain.values(), key=lambda r: r["delta_sasa"])
        record = dict(representative)
        record["per_chain"] = {
            cid: {
                "delta_sasa": r["delta_sasa"],
                "min_distance": r["min_distance"],
                "contacts": len(r["partners"]),
            }
            for cid, r in sorted(by_chain.items())
        }
        record["in_contact"] = record["min_distance"] is not None
        record["contacts_glycan_only"] = bool(record["partners"]) and all(
            p["is_glycan"] for p in record["partners"]
        )
        positions[str(seq_id)] = record

    glycan_records.sort(key=lambda r: (r["glycan_chain"], r["seq"]))
    return positions, glycan_records


def main() -> None:
    out = {
        "_schema": "fcatlas-interface/1",
        "_about": (
            "Per-position interface detail for the deposited complexes used by the "
            "atlas: buried surface area, named partner residues, and the "
            "provenance each deposition states about its own chains."
        ),
        "_method": {
            "sasa_algorithm": "Shrake-Rupley",
            "probe_radius_angstrom": PROBE,
            "test_points_per_atom": N_POINTS,
            "atoms": "heavy atoms only; waters and ions removed",
            "glycans": (
                "glycan-only chains are adopted by the nearest protein chain, so an "
                "Fc glycan is Fc surface for both area and contacts"
            ),
            "sasa_free": "area of the Fc chains alone, including their glycans",
            "sasa_bound": "area of the same residues in the complex",
            "delta_sasa": "sasa_free minus sasa_bound, square angstroms",
            "buried_fraction": "delta_sasa divided by sasa_free",
            "chain_choice": (
                "each position is reported for the Fc chain that buries the most "
                "surface, since receptor binding to the Fc dimer is asymmetric; "
                "every chain is kept under per_chain"
            ),
            "contact_cutoff_angstrom": CUTOFF,
            "min_delta_sasa_reported": MIN_DELTA,
            "numbering": "auth_seq_id, which is EU numbering in these entries",
        },
        "complexes": {},
    }

    for pdb, spec in JOBS.items():
        model = load_model(cif_path(pdb))
        meta = metadata(cif_path(pdb))
        fc_glycan, partner_glycan, ignored = classify_chains(
            model, spec["fc"], spec["partner"]
        )
        fc_all = spec["fc"] + fc_glycan
        partner_all = spec["partner"] + partner_glycan

        free = sasa_by_residue(model, fc_all)
        bound = sasa_by_residue(model, fc_all + partner_all)
        contacts = partner_contacts(model, fc_all, partner_all)
        positions, glycan_records = build_positions(
            free, bound, contacts, spec["fc"], fc_glycan
        )

        fc_entities = [
            e for e in meta["entities"] if set(e["chains"]) & set(spec["fc"])
        ]
        partner_entities = [
            e for e in meta["entities"] if set(e["chains"]) & set(spec["partner"])
        ]

        out["complexes"][pdb] = {
            "method": meta["method"],
            "resolution_angstrom": meta["resolution_angstrom"],
            "fc_chains": spec["fc"],
            "fc_glycan_chains": fc_glycan,
            "partner_chains": spec["partner"],
            "partner_glycan_chains": partner_glycan,
            "ignored_chains": ignored,
            "deposited_fc_mutation": next(
                (e["deposited_mutation"] for e in fc_entities), None
            ),
            "deposited_partner_mutation": next(
                (e["deposited_mutation"] for e in partner_entities), None
            ),
            "entities": meta["entities"],
            "positions": positions,
            "fc_glycan_contacts": glycan_records,
        }

        contact_n = sum(1 for r in positions.values() if r["in_contact"])
        buried_only = len(positions) - contact_n
        total = sum(r["delta_sasa"] for r in positions.values()) + sum(
            r["delta_sasa"] for r in glycan_records
        )
        print(
            f"{pdb}: {contact_n} positions in contact, {buried_only} buried without "
            f"contact, {len(glycan_records)} Fc glycan residues involved, "
            f"{total:.0f} A^2 buried; fc mutation="
            f"{out['complexes'][pdb]['deposited_fc_mutation']!r} partner mutation="
            f"{out['complexes'][pdb]['deposited_partner_mutation']!r}"
        )

    with open(REPO / "data" / "interface_detail.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print("wrote interface_detail.json")


if __name__ == "__main__":
    main()
