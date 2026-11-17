"""Izris grafov za poročilo iz rezultatov analyze.py.

Uporaba:
    python plot.py --summary runs/<izpad>/summary.csv runs/<premik>/summary.csv --out img

`summary.csv` je lahko en sam (zagon `all`) ali več (ločena zagona `removal` in `shift`).
Izriše dva grafa v vektorskem PDF:
  okvare_postavitev.pdf  izpad in premik posamezne celice, prikazano na postavitvi vrat
  premik_razdalja.pdf    delež premikov, pri katerih vrata odpovejo, glede na razdaljo
Postavitev vrat in valovne oblike kontrole so posnetki zaslona iz QCADesigner-E.
Zahteva matplotlib.
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parent
CM = 1 / 2.54
TEXT_W = 11.7 * CM  # širina besedila v svmult
CELL = 18  # širina celice [nm]

INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
ZONE = {"0": "#b7d3f6", "1": "#f9c9ad"}  # urina cona 0 (modra), 1 (oranžna)
ZONE_EDGE = {"0": "#2a78d6", "1": "#eb6834"}
OK_FILL, FAIL_FILL, NA_FILL = "#cde2fb", "#eb6834", "#efeeeb"
SEQ = LinearSegmentedColormap.from_list(
    "seq", ["#f4f8fd", "#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def load_cells():
    with open(ROOT / "cells.csv", newline="") as f:
        cells = list(csv.DictReader(f))
    for c in cells:
        c["x"], c["y"] = float(c["x_nm"]), float(c["y_nm"])
    return cells


def load_summary(paths):
    rows = {}
    for p in paths:
        with open(p, newline="") as f:
            for r in csv.DictReader(f):
                rows[r["variant"]] = r
    return rows


def n_wrong(row):
    return len(row["wrong_combos_ABCDE"].split())


def draw_layout(ax, cells, fill, edge=None, lw=None, small=None, tcolor=None):
    """Celice narišemo na njihovih položajih (os y navzdol), z oznako iz članka in dodatnim napisom."""
    for c in cells:
        lab = c["paper_label"]
        ax.add_patch(Rectangle((c["x"] - CELL / 2, c["y"] - CELL / 2), CELL, CELL, facecolor=fill(c),
                               linewidth=(lw or {}).get(lab, 0.6), edgecolor=(edge or {}).get(lab, INK2)))
        t = {"OUT": "izhod", "E": "F"}.get(lab, lab)  # vhod E iz članka je v MV5 in GUI poimenovan F
        col = (tcolor or {}).get(lab, INK)
        ax.text(c["x"], c["y"] - (2.6 if small and lab in small else 0), t, ha="center", va="center",
                fontsize=6.5 if len(t) <= 2 else 5.5, fontweight="normal" if c["role"] == "normal" else "bold",
                color=col)
        if small and lab in small:
            ax.text(c["x"], c["y"] + 4.4, small[lab], ha="center", va="center", fontsize=5.5, color=col)
    xs = [c["x"] for c in cells]
    ys = [c["y"] for c in cells]
    ax.set_xlim(min(xs) - 12, max(xs) + 12)
    ax.set_ylim(max(ys) + 12, min(ys) - 12)
    ax.set_aspect("equal")
    ax.axis("off")


def fig_faults_layout(cells, summary, out):
    """(a) izid izpada vsake notranje celice, (b) število premikov celice (od 12), pri katerih vrata odpovejo."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(TEXT_W, 7.4 * CM))

    # (a) izpad; članek (Tabela 4): vrata odpovejo pri izpadu celic 01, 02, 04 in 08
    paper_fail = {"01", "02", "04", "08"}
    state, small, edge, lw = {}, {}, {}, {}
    for c in cells:
        lab = c["paper_label"]
        r = summary.get(f"remove_{lab}")
        if r is None:
            state[lab] = "na"
        elif r["verdict"] == "OK":
            state[lab] = "ok"
        else:
            state[lab] = "fail"
            small[lab] = f"{n_wrong(r)}/32"
            if lab not in paper_fail:  # razhajanje s člankom
                edge[lab], lw[lab] = INK, 1.8
    fills = {"ok": OK_FILL, "fail": FAIL_FILL, "na": NA_FILL}
    draw_layout(a, cells, lambda c: fills[state[c["paper_label"]]], edge=edge, lw=lw, small=small)
    a.set_title("(a) izpad celice", fontsize=8, pad=3)
    a.legend(handles=[Patch(facecolor=OK_FILL, edgecolor=INK2, linewidth=0.6, label="vrata delujejo"),
                      Patch(facecolor=FAIL_FILL, edgecolor=INK2, linewidth=0.6, label="vrata odpovejo"),
                      Patch(facecolor=FAIL_FILL, edgecolor=INK, linewidth=1.8,
                            label="odpovejo, v članku delujejo"),
                      Patch(facecolor=NA_FILL, edgecolor=INK2, linewidth=0.6, label="ni preizkušeno")],
             loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=1, frameon=False, handlelength=1.2,
             borderaxespad=0)

    # (b) premik: 4 smeri x 3 razdalje = 12 premikov na celico
    fails = {c["paper_label"]: 0 for c in cells}
    for r in summary.values():
        if r["type"] == "shift" and r["verdict"] != "OK":
            fails[r["paper_label"]] += 1
    vmax = 12
    draw_layout(b, cells, lambda c: SEQ(0.03 + 0.97 * fails[c["paper_label"]] / vmax) if fails[c["paper_label"]]
                else "#ffffff", small={k: f"{v}/12" for k, v in fails.items()},
                tcolor={k: "#ffffff" for k, v in fails.items() if v >= 7})
    b.set_title("(b) premik celice", fontsize=8, pad=3)
    sm = plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, vmax))
    cax = b.inset_axes([0.1, -0.13, 0.8, 0.035])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[0, 3, 6, 9, 12])
    cb.outline.set_linewidth(0.4)
    cb.set_label("premiki, pri katerih vrata odpovejo", fontsize=7, labelpad=2)
    fig.subplots_adjust(wspace=0.05)
    fig.savefig(out / "okvare_postavitev.pdf")
    plt.close(fig)


def fig_distance(cells, summary, out):
    """Delež premikov, pri katerih vrata odpovejo, glede na razdaljo in del vezja."""
    zone = {c["paper_label"]: c["clock"] for c in cells}
    stats = {}
    for r in summary.values():
        if r["type"] != "shift":
            continue
        s = stats.setdefault((float(r["distance_nm"]), zone[r["paper_label"]]), {"n": 0, "f": 0, "fov": 0})
        s["n"] += 1
        if r["verdict"] != "OK":
            s["f"] += 1
            s["fov"] += r["overlap"] == "1"
    dists = sorted({d for d, _ in stats})
    fig, ax = plt.subplots(figsize=(8.0 * CM, 5.0 * CM))
    w = 0.36
    for k, z in enumerate(("0", "1")):
        for i, d in enumerate(dists):
            s = stats[(d, z)]
            x = i + (k - 0.5) * (w + 0.04)
            clean = 100 * (s["f"] - s["fov"]) / s["n"]
            over = 100 * s["fov"] / s["n"]
            ax.bar(x, clean, w, color=ZONE_EDGE[z], linewidth=0)
            ax.bar(x, over, w, bottom=clean, color=ZONE[z], hatch="////", edgecolor=ZONE_EDGE[z], linewidth=0)
            ax.text(x, clean + over + 2, f"{s['f']}/{s['n']}", ha="center", va="bottom", fontsize=6.5, color=INK2)
    ax.set_xticks(range(len(dists)))
    ax.set_xticklabels([f"{d:g} nm" for d in dists])
    ax.set_xlabel("razdalja premika")
    ax.set_ylabel("vrata odpovejo [%]")
    ax.set_ylim(0, 108)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.yaxis.grid(True, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.legend(handles=[Patch(color=ZONE_EDGE["0"], label="vhodni del (cona 0)"),
                       Patch(color=ZONE_EDGE["1"], label="izhodni del (cona 1)"),
                       Patch(facecolor="#ffffff", edgecolor=INK2, hatch="////", linewidth=0,
                             label="od tega s prekrivanjem")],
              loc="upper left", frameon=False, handlelength=1.2)
    fig.savefig(out / "premik_razdalja.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", nargs="+", required=True, help="summary.csv iz analyze.py")
    ap.add_argument("--out", default="img")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cells = load_cells()
    summary = load_summary(args.summary)
    fig_faults_layout(cells, summary, out)
    fig_distance(cells, summary, out)
    print(f"grafa zapisana v {out}")


if __name__ == "__main__":
    main()
