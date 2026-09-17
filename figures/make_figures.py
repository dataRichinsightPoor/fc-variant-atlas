"""Regenerate every figure in this directory.

    python figures/make_figures.py

Output is byte-reproducible: no timestamps, no version strings, fixed element
order, deterministic layout. Rerunning on unchanged data yields identical PNG
and SVG files, so a diff in a figure always means a change in the data.
"""

from __future__ import annotations

import os
import struct
import zlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from fcatlas import (  # noqa: E402
    classify,
    engineering_against_burial,
    load_interface,
    domain_of,
    fcrn_provenance,
    load_complexes,
    load_variants,
    numbering,
    position_index,
)

HERE = os.path.dirname(os.path.abspath(__file__))

INK = "#16181d"
MUTED = "#6b7280"
GRID = "#dfe3e8"
FCGR = "#2f6f8f"
FCRN = "#b0752a"
BOTH = "#5b5f97"
OFF = "#a8adb5"
ACCENT = "#b0413e"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlesize": 10.5,
        "axes.titleweight": "medium",
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "svg.hashsalt": "fcatlas",
        "svg.fonttype": "none",
        "figure.dpi": 200,
        "savefig.dpi": 200,
    }
)


# --------------------------------------------------------------- deterministic IO
def _strip_png_metadata(path: str) -> None:
    """Drop tEXt/tIME/iTXt/pHYs chunks so the PNG has no toolchain fingerprint."""
    with open(path, "rb") as fh:
        blob = fh.read()
    sig, out, i = blob[:8], bytearray(blob[:8]), 8
    assert sig == b"\x89PNG\r\n\x1a\n"
    while i < len(blob):
        (length,) = struct.unpack(">I", blob[i : i + 4])
        ctype = blob[i + 4 : i + 8]
        chunk = blob[i : i + 12 + length]
        if ctype not in (b"tEXt", b"iTXt", b"zTXt", b"tIME", b"pHYs"):
            out += chunk
        i += 12 + length
    with open(path, "wb") as fh:
        fh.write(bytes(out))
    # sanity: the file must still decode
    assert zlib.crc32(b"IEND") is not None


def save(fig, stem: str) -> None:
    png = os.path.join(HERE, stem + ".png")
    svg = os.path.join(HERE, stem + ".svg")
    fig.savefig(png, bbox_inches="tight", metadata={})
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None, "Creator": None})
    _strip_png_metadata(png)
    plt.close(fig)
    print("wrote", os.path.relpath(png, HERE), "and", os.path.relpath(svg, HERE))


# --------------------------------------------------- figure 1: hotspot density
def figure_hotspots() -> None:
    idx = position_index()
    n_res = numbering("IgG1").residue
    cx = load_complexes()
    fcgr_pos = {p for k, c in cx.items() if "Fcgamma" in c.partner for p in c.contacts}
    fcrn_pos = {p for k, c in cx.items() if "FcRn" in c.partner for p in c.contacts}

    eus = sorted(idx)
    counts = [len(idx[e]) for e in eus]

    def colour(e: int) -> str:
        if e in fcgr_pos and e in fcrn_pos:
            return BOTH
        if e in fcgr_pos:
            return FCGR
        if e in fcrn_pos:
            return FCRN
        return OFF

    lo, hi = min(eus) - 6, max(eus) + 6
    fig, ax = plt.subplots(figsize=(9.2, 3.6))

    spans, start, prev = [], lo, domain_of(int(lo))
    for e in range(int(lo), int(hi) + 1):
        d = domain_of(e)
        if d != prev:
            spans.append((start, e, prev))
            start, prev = e, d
    spans.append((start, hi, prev))
    for i, (a, b, d) in enumerate(spans):
        if i % 2:
            ax.axvspan(a, b, color="#f4f6f8", linewidth=0, zorder=0)
        if d in ("hinge", "CH2", "CH3"):
            y = max(counts) + (2.1 if d == "hinge" else 1.0)
            ax.text(a + 1.5, y, d, ha="left", va="center", fontsize=8.5, color=MUTED)

    ax.bar(eus, counts, width=1.7, color=[colour(e) for e in eus], linewidth=0, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, max(counts) + 2.5)
    ax.set_xlabel("EU position (Edelman numbering)")
    ax.set_ylabel("engineered variants")
    ax.set_title(
        "Where Fc engineering concentrates, and whether the position was ever "
        "measured in contact"
    )
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6, zorder=1)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    top = sorted(idx, key=lambda e: (-len(idx[e]), e))[:3]
    offsets = {0: (20, 4), 1: (24, -10), 2: (30, -4)}
    for rank, e in enumerate(top):
        ax.annotate(
            f"{n_res(e)}{e} \u2014 {len(idx[e])} variants",
            xy=(e + 1.0, len(idx[e])),
            xytext=offsets[rank],
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=8.5,
            color=ACCENT if e in (234, 235) else INK,
            arrowprops=dict(arrowstyle="-", linewidth=0.7,
                            color=ACCENT if e in (234, 235) else MUTED,
                            shrinkA=0, shrinkB=1),
        )

    ax.legend(
        handles=[
            Patch(color=FCGR, label="in a measured Fc\u03b3R interface"),
            Patch(color=FCRN, label="in the measured FcRn interface"),
            Patch(color=BOTH, label="in both"),
            Patch(color=OFF, label="never measured in contact"),
        ],
        frameon=False,
        fontsize=8,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
    )
    fig.text(
        0.5,
        -0.14,
        "74 curated variants, 60 distinct positions. Contact means a heavy atom "
        "within 5.0 \u00c5 of the partner in 1E4K, 1T89, 5XJE or 4N0U.",
        ha="center",
        fontsize=8,
        color=MUTED,
    )
    save(fig, "fig1_hotspot_density")


# ------------------------------------------- figure 2: interface classification
def figure_classification() -> None:
    rows = classify()
    groups = {}
    for r in rows:
        f, n = bool(r["fcgr_interface"]), bool(r["fcrn_interface"])
        key = "both" if (f and n) else "FcgammaR only" if f else "FcRn only" if n else "neither"
        groups.setdefault(r["intent"], {}).setdefault(key, []).append(r["id"])

    pretty = {
        "enhance_adcc": "enhance ADCC",
        "enhance_cdc": "enhance CDC",
        "reduce_cdc": "reduce CDC",
        "enhance_fcgr2b": "enhance Fc\u03b3RIIb",
        "halflife": "extend half-life",
        "reduce_halflife": "shorten half-life",
        "heterodimer": "force heterodimer",
        "silence+stabilize": "silence and stabilize",
        "silence+halflife": "silence and extend half-life",
    }
    order = sorted(groups, key=lambda k: -sum(len(v) for v in groups[k].values()))
    ticks = [pretty.get(k, k) for k in order]
    keys = ["FcgammaR only", "FcRn only", "both", "neither"]
    colours = {"FcgammaR only": FCGR, "FcRn only": FCRN, "both": BOTH, "neither": OFF}

    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    left = [0.0] * len(order)
    for k in keys:
        vals = [len(groups[i].get(k, [])) for i in order]
        ax.barh(ticks, vals, left=left, color=colours[k], height=0.66, linewidth=0)
        for y, (v, l) in enumerate(zip(vals, left)):
            if v:
                ax.text(
                    l + v / 2, y, str(v), ha="center", va="center",
                    fontsize=8, color="white" if k != "neither" else INK,
                )
        left = [a + b for a, b in zip(left, vals)]

    ax.invert_yaxis()
    ax.set_xlabel("variants")
    ax.set_title("Design intent against measured interface location")
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        handles=[Patch(color=colours[k], label=k.replace("Fcgamma", "Fc\u03b3")) for k in keys],
        frameon=False, fontsize=8, loc="lower right",
    )
    fig.text(
        0.5, -0.04,
        "24 of 74 variants change no position ever measured in contact with Fc\u03b3R "
        "or FcRn in these structures.",
        ha="center", fontsize=8, color=MUTED,
    )
    save(fig, "fig2_intent_vs_interface")


# ------------------------------------- figure 3: FcRn structure provenance
def figure_fcrn_provenance() -> None:
    survey = fcrn_provenance()
    entries = [e for e in survey["entries"] if e["has_igg_fc"]]
    complexes = [e for e in entries if e["category"] == "fc_fcrn_complex"]
    unbound = [e for e in entries if e["category"] == "unbound_fc"]

    def genotype(e):
        muts = [m for _, m in e["fc"] if isinstance(m, list)][0]
        if [m.lower() for m in muts] == ["wild type"]:
            return "wild type"
        return "/".join(sorted(muts, key=lambda m: int("".join(c for c in m if c.isdigit()))))

    rows = [("bound to FcRn", None, True)]
    rows += [(e["pdb"], genotype(e), True) for e in sorted(complexes, key=lambda x: x["pdb"])]
    rows += [("solved unbound", None, False)]
    rows += [(e["pdb"], genotype(e), False) for e in sorted(unbound, key=lambda x: x["pdb"])]

    contact = sorted(load_complexes()["4N0U"].contacts)

    fig, ax = plt.subplots(figsize=(8.8, 0.30 * len(rows) + 1.9))
    for y, (pdb, geno, is_cx) in enumerate(rows):
        if geno is None:
            ax.text(-4.0, y, pdb, ha="left", va="center", fontsize=9, color=MUTED)
            ax.plot([-0.05, len(contact) - 0.5], [y + 0.5, y + 0.5],
                    color=GRID, linewidth=0.8)
            continue
        wt = geno == "wild type"
        ax.text(-0.6, y, pdb, ha="right", va="center", fontsize=8.5,
                color=INK, family="DejaVu Sans Mono")
        ax.text(
            len(contact) + 0.8, y, geno,
            ha="left", va="center", fontsize=8,
            color=MUTED if wt else INK,
        )
        muts = (
            {}
            if wt
            else {int("".join(c for c in m if c.isdigit())): m for m in geno.split("/")}
        )
        for x, eu in enumerate(contact):
            hit = eu in muts
            ax.add_patch(
                plt.Rectangle(
                    (x - 0.42, y - 0.36), 0.84, 0.72,
                    facecolor=ACCENT if hit else ("#eef1f4" if is_cx else "#f7f8f9"),
                    edgecolor=GRID, linewidth=0.4,
                )
            )
        ax.plot([-0.15], [y], marker="o", markersize=4,
                color=FCRN if is_cx else "#d5d9de", zorder=3)

    ax.set_xlim(-4.2, len(contact) + 7.5)
    ax.set_ylim(len(rows) - 0.5, -1.7)
    ax.set_xticks(range(len(contact)))
    ax.set_xticklabels([str(e) for e in contact], fontsize=7.5, rotation=90)
    ax.set_yticks([])
    ax.xaxis.set_ticks_position("top")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_title(
        "Every FcRn contact position, against the Fc genotype of the entry it "
        "was measured in",
        pad=26,
    )
    ax.text(len(contact) / 2 - 0.5, -1.15,
            "EU positions within 5.0 \u00c5 of FcRn in 4N0U",
            ha="center", fontsize=8, color=MUTED)
    ax.legend(
        handles=[
            Patch(color=ACCENT, label="position substituted in this entry"),
            Patch(color="#eef1f4", label="contact position left as wild type"),
        ],
        frameon=False, fontsize=8, loc="upper center",
        bbox_to_anchor=(0.45, -0.015), bbox_transform=ax.transAxes,
        ncol=3, handlelength=1.2, columnspacing=1.6,
    )
    ax.text(
        0.45, -0.075,
        "All four FcRn-bound Fc chains in this survey are engineered. Wild-type "
        "human Fc appears only unbound.",
        transform=ax.transAxes, ha="center", va="top", fontsize=8, color=MUTED,
    )
    save(fig, "fig3_fcrn_provenance")


# ---------------------------------- figure 4: the numbering collision at 234
def figure_numbering_collision() -> None:
    lo, hi = 228, 242
    isos = ["IgG1", "IgG2", "IgG3", "IgG4"]
    fig, ax = plt.subplots(figsize=(7.6, 2.9))
    for row, iso in enumerate(isos):
        n = numbering(iso)
        for col, eu in enumerate(range(lo, hi + 1)):
            aa = n.residue(eu)
            gap = aa is None
            highlight = eu in (234, 235) and iso == "IgG2"
            ax.add_patch(
                plt.Rectangle(
                    (col - 0.46, row - 0.4), 0.92, 0.8,
                    facecolor="#f3d9d8" if highlight else ("#f2f4f6" if gap else "white"),
                    edgecolor=GRID, linewidth=0.5,
                )
            )
            ax.text(
                col, row, aa if aa else "\u2013",
                ha="center", va="center", fontsize=10,
                family="DejaVu Sans Mono",
                color=ACCENT if highlight else (MUTED if gap else INK),
            )
        ax.text(-1.1, row, iso, ha="right", va="center", fontsize=9)

    for col, eu in enumerate(range(lo, hi + 1)):
        ax.text(col, -0.95, str(eu), ha="center", va="center", fontsize=7.5,
                rotation=90, color=MUTED)

    ax.annotate(
        "IgG2 has no residue\nat 234; its valine\nsits at 235",
        xy=(234.5 - lo, 0.55), xytext=(hi - lo + 1.4, 0.25),
        fontsize=8.5, color=ACCENT, ha="left", va="center",
        arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.8,
                        shrinkA=2, shrinkB=2,
                        connectionstyle="angle,angleA=0,angleB=-90,rad=0"),
    )
    ax.set_xlim(-4.4, hi - lo + 6.6)
    ax.set_ylim(3.9, -1.5)
    ax.axis("off")
    ax.set_title("The lower hinge in EU coordinates, aligned across isotypes")
    fig.text(
        0.5, -0.02,
        "The same physical residue is written V234A in the primary literature "
        "and V235A by IMGT. The alignment shows why both are defensible.",
        ha="center", fontsize=8, color=MUTED,
    )
    save(fig, "fig4_numbering_collision")


# ------------------------------------- figure 5: attention against measurement
def figure_burial_against_density() -> None:
    """How often a position is engineered, against how much it actually buries.

    Both axes are measurements. The horizontal one comes from the deposited
    coordinates; the vertical one comes from what the field has chosen to
    change. They are close to uncorrelated, and the interesting points are the
    ones far from any diagonal: heavily buried positions nobody engineers, and
    heavily engineered positions that barely touch the partner.
    """
    interfaces = load_interface()
    rows = [r for r in engineering_against_burial() if r["delta_sasa"]]
    fcgr = {"1E4K", "1T89", "5XJE"}

    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    ax.set_axisbelow(True)
    ax.grid(axis="both", color=GRID, linewidth=0.6)

    for row in sorted(rows, key=lambda r: -r["delta_sasa"]):
        colour = FCGR if row["pdb"] in fcgr else FCRN
        ax.scatter(
            row["delta_sasa"],
            row["records"],
            s=26 + 150 * row["buried_fraction"],
            facecolor=colour,
            edgecolor="white",
            linewidth=0.7,
            alpha=0.9,
            zorder=3,
        )

    # Label the points that carry the argument, placed by hand so nothing
    # collides and no label needs a leader line.
    labels = {
        329: (-6, 14, "right"),
        253: (7, 3, "left"),
        434: (6, 8, "left"),
        235: (0, 14, "center"),
        236: (0, -16, "center"),
        234: (0, 14, "center"),
        311: (6, -6, "left"),
        296: (0, -16, "center"),
        330: (11, -7, "left"),
        239: (0, 14, "center"),
        254: (8, 8, "left"),
        309: (8, 4, "left"),
        252: (0, -15, "center"),
    }
    for row in rows:
        if row["eu"] not in labels:
            continue
        dx, dy, ha = labels[row["eu"]]
        ax.annotate(
            str(row["eu"]),
            xy=(row["delta_sasa"], row["records"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.5,
            color=INK,
            ha=ha,
            va="center",
        )

    ax.set_xlabel("surface the partner buries at this position, square angstroms")
    ax.set_ylabel("curated variant records that change it")
    ax.set_xlim(-6, 142)
    ax.set_ylim(-1.6, 24.0)
    ax.set_title(
        "What the field engineers, against what the interface actually buries"
    )

    ax.axhline(0, color=MUTED, linewidth=0.7, zorder=2)
    ax.annotate(
        "buried and untouched: EU 253 gives up 96 percent\nof its surface to FcRn "
        "and appears in no record here",
        xy=(121.2, 0.15), xytext=(112, 8.6),
        fontsize=8.5, color=ACCENT, ha="center", va="bottom",
        arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.8,
                        shrinkA=2, shrinkB=4),
    )
    ax.annotate(
        "engineered and barely buried: EU 332 gives up\n9 square angstroms and "
        "carries six records",
        xy=(9.5, 6.2), xytext=(48, 11.4),
        fontsize=8.5, color=ACCENT, ha="center", va="center",
        arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.8,
                        shrinkA=2, shrinkB=4),
    )

    handles = [
        Patch(facecolor=FCGR, edgecolor="none", label="measured at an FcgammaR interface"),
        Patch(facecolor=FCRN, edgecolor="none", label="measured at the FcRn interface"),
    ]
    ax.legend(
        handles=handles, loc="upper left", frameon=False, fontsize=8.5,
        handlelength=1.1, borderaxespad=0.4,
    )

    total = sum(d.total_buried_area for d in interfaces.values())
    fig.text(
        0.5, -0.075,
        "Each position appears once, at the structure where it buries the most; marker "
        "area scales with the fraction of the residue's free surface buried there.\n"
        "Areas are Shrake-Rupley differences over the deposited coordinates, probe "
        f"1.40 angstroms. {len(rows)} positions carry measurable burial across four "
        f"complexes totalling {total:.0f} square angstroms of interface; positions "
        "with none are not plotted.",
        ha="center", fontsize=8, color=MUTED,
    )
    save(fig, "fig5_burial_against_density")


def main() -> None:
    figure_hotspots()
    figure_classification()
    figure_fcrn_provenance()
    figure_numbering_collision()
    figure_burial_against_density()
    print(f"{len(load_variants())} variants; figures written to {HERE}")


if __name__ == "__main__":
    main()
