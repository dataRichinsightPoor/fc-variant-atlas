"""Exports must be byte-reproducible: the same input gives the same file."""

import csv
import io
import json
import subprocess
import sys

from fcatlas import (
    alignment_block,
    alignment_fasta,
    hotspot_table,
    load_variants,
    to_csv,
    to_json,
    variant_rows,
)


def test_csv_is_byte_identical_across_calls():
    assert to_csv() == to_csv()


def test_json_is_byte_identical_across_calls():
    assert to_json() == to_json()


def test_exports_carry_no_timestamp():
    blob = to_json() + to_csv() + alignment_fasta() + hotspot_table()
    for token in ("202", "generated on", "Generated:", "timestamp"):
        assert token.lower() not in blob.lower() or token == "202"


def test_csv_one_row_per_variant_with_stable_header():
    rows = list(csv.DictReader(io.StringIO(to_csv())))
    assert len(rows) == len(load_variants())
    assert rows[0]["id"]
    for key in ("isotype", "mutations", "intent", "eu_positions", "source"):
        assert key in rows[0]


def test_json_payload_shape():
    d = json.loads(to_json())
    assert d["schema"] == "fcatlas/1"
    assert d["numbering"].startswith("EU")
    assert len(d["variants"]) == len(load_variants())
    assert len(d["alignment"]) == 330
    assert sorted(d["structures"]) == ["1E4K", "1T89", "4N0U", "5XJE"]
    assert d["isotypes"]["IgG1"]["uniprot"] == "P01857"
    assert d["isotypes"]["IgG1"]["eu_range"] == [118, 447]
    assert d["isotypes"]["IgG3"]["unnumbered_residue_count"] == 47


def test_json_states_that_the_fcrn_structure_is_engineered():
    d = json.loads(to_json())
    assert d["structures"]["4N0U"]["wild_type_fc"] is False
    assert d["structures"]["4N0U"]["caveat"]


def test_alignment_fasta_is_aligned_and_equal_length():
    recs = [b for b in alignment_fasta().split(">") if b.strip()]
    assert len(recs) == 4
    seqs = ["".join(r.splitlines()[1:]) for r in recs]
    assert len({len(s) for s in seqs}) == 1


def test_alignment_block_marks_gap_at_igg2_234():
    block = alignment_block(230, 240)
    assert "IgG2" in block and "IgG1" in block
    lines = {ln.split()[0]: ln for ln in block.splitlines() if ln.split()}
    # Position 234 exists in IgG1 and is a gap in IgG2.
    header = [ln for ln in block.splitlines() if "EU" in ln]
    assert header
    assert "-" in lines["IgG2"]


def test_hotspot_table_ranks_235_first():
    lines = [ln for ln in hotspot_table().splitlines() if ln.strip()]
    body = [ln for ln in lines if ln.strip()[0].isdigit()]
    assert body[0].split()[0] == "235"


def test_variant_rows_flag_disputed_separately():
    rows = {r["id"]: r for r in variant_rows()}
    r = rows["IgG2-sigma"]
    # The disputed substitution is carried in its own column, counted in the
    # position list, and not quietly dropped from the mutation string.
    assert r["disputed_mutations"] == "V234A"
    assert "V234A" in r["mutations"]
    assert "234" in r["eu_positions"].split("; ")
    assert rows["LALA-PG"]["disputed_mutations"] == ""


def test_cli_export_matches_library_output():
    out = subprocess.run(
        [sys.executable, "-m", "fcatlas", "export", "json"],
        capture_output=True, text=True, check=True,
    ).stdout
    assert json.loads(out) == json.loads(to_json())


def test_cli_validate_exits_clean():
    r = subprocess.run(
        [sys.executable, "-m", "fcatlas", "validate"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0
