"""Fault sweeps for the MV5 QCA gate (see PLAN.md).

Generates design variants (single-cell removal, single-cell shift), simulates
each one headless with engine/bin/qca_sim inside WSL, records everything
needed to repeat the run in manifest.json, then runs analyze.py on it.

    python sweep.py removal
    python sweep.py shift --distances 1 2 5
    python sweep.py all --seed 2 --name seed2

Each run goes to runs/<date>_<time>_<kind>[_<name>]/ unless --out is given.
"""
import argparse
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import analyze

ROOT = Path(__file__).resolve().parent
QCA_SIM = ROOT / "engine" / "bin" / "qca_sim"
BUILD_SH = ROOT / "engine" / "build.sh"
CELLS_CSV = ROOT / "cells.csv"
CELL_SIZE = 18.0

# Screen coordinates: y grows downwards
DIRECTIONS = {"up": (0, -1), "down": (0, 1), "left": (-1, 0), "right": (1, 0)}

CELL_BLOCK = re.compile(r"\[TYPE:QCADCell\].*?\[#TYPE:QCADCell\]\n", re.S)
# Absolute coordinates inside a cell block (centre, bounding box, dots, label);
# cxWorld/cyWorld are sizes and must not move.
ABS_COORD = re.compile(r"^(x|y|bounding_box\.xWorld|bounding_box\.yWorld)=(-?[0-9.]+)$", re.M)


def to_wsl(path):
    p = Path(path).resolve()
    return "/mnt/" + p.drive[0].lower() + p.as_posix()[2:]


def wsl(script):
    return subprocess.run(["wsl", "-d", "Ubuntu", "--", "bash", "-c", script], capture_output=True, text=True)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalise_decimals(text):
    # Files saved under a comma-decimal locale ("440,000000"); the loader
    # (g_ascii_strtod) stops at the comma.
    return re.sub(r"^([A-Za-z_.]+=-?[0-9]+),([0-9]+)", r"\1.\2", text, flags=re.M)


def load_cells():
    import csv
    with open(CELLS_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["file_index"] = int(r["file_index"])
        r["x_nm"] = float(r["x_nm"])
        r["y_nm"] = float(r["y_nm"])
        r["clock"] = int(r["clock"])
    return sorted(rows, key=lambda r: r["file_index"])


def split_design(text):
    blocks = list(CELL_BLOCK.finditer(text))
    return text[: blocks[0].start()], [b.group(0) for b in blocks], text[blocks[-1].end():]


def join_design(head, blocks, tail):
    return head + "".join(blocks) + tail


def shift_block(block, dx, dy):
    def repl(m):
        key, val = m.group(1), float(m.group(2))
        val += dx if key in ("x", "bounding_box.xWorld") else dy
        return f"{key}={val:.6f}"
    return ABS_COORD.sub(repl, block)


def min_gap(cells, moved_index, dx, dy):
    """Smallest edge-to-edge gap (nm) between the moved cell and any other; <0 means overlap."""
    me = cells[moved_index]
    mx, my = me["x_nm"] + dx, me["y_nm"] + dy
    gaps = [max(abs(c["x_nm"] - mx) - CELL_SIZE, abs(c["y_nm"] - my) - CELL_SIZE)  # apart along either axis
            for c in cells if c["file_index"] != moved_index]
    return min(gaps)


def build_variants(kind, design_text, cells, distances):
    """Returns [(metadata, design_text)] in a fixed order."""
    head, blocks, tail = split_design(design_text)
    if len(blocks) != len(cells):
        sys.exit(f"{len(blocks)} cell blocks in the design vs {len(cells)} rows in cells.csv")
    variants = [({"name": "control", "type": "control"}, design_text)]

    def cell_meta(c):
        return {"paper_label": c["paper_label"], "file_index": c["file_index"], "role": c["role"]}

    if kind in ("removal", "all"):
        for c in cells:
            if c["role"] != "normal":
                continue
            kept = [b for i, b in enumerate(blocks) if i != c["file_index"]]
            meta = {"name": f"remove_{c['paper_label']}", "type": "removal", **cell_meta(c)}
            variants.append((meta, join_design(head, kept, tail)))

    if kind in ("shift", "all"):
        for c in cells:
            for dname, (ux, uy) in DIRECTIONS.items():
                for d in distances:
                    dx, dy = ux * d, uy * d
                    moved = list(blocks)
                    moved[c["file_index"]] = shift_block(blocks[c["file_index"]], dx, dy)
                    gap = min_gap(cells, c["file_index"], dx, dy)
                    meta = {"name": f"shift_{c['paper_label']}_{dname}_{d:g}nm", "type": "shift", **cell_meta(c),
                            "direction": dname, "distance_nm": d, "min_gap_nm": gap, "overlap": gap < 0}
                    variants.append((meta, join_design(head, moved, tail)))
    return variants


def ensure_binary(rebuild):
    if rebuild or not QCA_SIM.exists():
        proc = wsl(f"bash {to_wsl(BUILD_SH)}")
        if proc.returncode != 0:
            sys.exit(f"build failed:\n{proc.stdout}\n{proc.stderr}")
        print(proc.stdout.strip().splitlines()[-1])


def simulate(run_dir, samples, seed, jobs):
    # One WSL call; xargs runs the simulations in parallel
    script = (
        f"cd {to_wsl(run_dir)} && ls designs | sed 's/\\.qca$//' | "
        f"xargs -P {jobs} -I{{}} sh -c '{to_wsl(QCA_SIM)} -f designs/{{}}.qca -n {samples} -s {seed} "
        f"-o raw/{{}}.csv >/dev/null 2>raw/{{}}.log || echo FAILED {{}}'"
    )
    proc = wsl(script)
    failed = sorted(l.split()[1] for l in proc.stdout.splitlines() if l.startswith("FAILED"))
    if proc.returncode != 0 or failed:
        sys.exit(f"simulation failed for: {failed or proc.stderr}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["removal", "shift", "all"])
    ap.add_argument("--design", default=str(ROOT / "MV5"))
    ap.add_argument("--distances", type=float, nargs="+", default=[1, 2, 5])
    ap.add_argument("--samples", type=int, default=128000, help="paper Table 1: 128000")
    ap.add_argument("--seed", type=int, default=1, help="seed for the engine's cell-order randomisation")
    ap.add_argument("--jobs", type=int, default=8, help="parallel simulations (does not affect results)")
    ap.add_argument("--unreadable", type=float, default=0.1, help="|P| below this counts as a wrong output")
    ap.add_argument("--weak-abs", type=float, default=0.5, help="WEAK if min|P| below this")
    ap.add_argument("--weak-rel", type=float, default=0.5, help="WEAK if min|P| below this fraction of control")
    ap.add_argument("--name", default="", help="suffix for the run folder name")
    ap.add_argument("--out", default=None, help="explicit run folder (overwritten)")
    ap.add_argument("--rebuild", action="store_true", help="rebuild qca_sim before running")
    args = ap.parse_args()

    now = datetime.now()
    run_name = f"{now:%Y-%m-%d_%H%M%S}_{args.kind}" + (f"_{args.name}" if args.name else "")
    run_dir = Path(args.out) if args.out else ROOT / "runs" / run_name
    for sub in ("designs", "raw"):
        shutil.rmtree(run_dir / sub, ignore_errors=True)
        (run_dir / sub).mkdir(parents=True)

    ensure_binary(args.rebuild)

    cells = load_cells()
    design_path = Path(args.design).resolve()
    text = normalise_decimals(design_path.read_text(encoding="utf-8", errors="replace"))
    variants = build_variants(args.kind, text, cells, args.distances)
    for meta, vtext in variants:
        (run_dir / "designs" / f"{meta['name']}.qca").write_text(vtext, encoding="utf-8", newline="\n")

    manifest = {
        "run_name": run_dir.name,
        "kind": args.kind,
        "created": now.isoformat(timespec="seconds"),
        "command": " ".join(["python"] + sys.argv),
        "design": design_path.name,
        "hashes": {"design": sha256(design_path), "qca_sim": sha256(QCA_SIM), "cells_csv": sha256(CELLS_CSV),
                   "sweep_py": sha256(__file__), "analyze_py": sha256(analyze.__file__)},
        "python": platform.python_version(),
        "params": {"engine": "BISTABLE", "sim_type": "EXHAUSTIVE_VERIFICATION", "samples": args.samples,
                   "seed": args.seed, "distances": args.distances},
        "thresholds": {"unreadable": args.unreadable, "weak_abs": args.weak_abs, "weak_rel": args.weak_rel},
        "cells": cells,
        "variants": [meta for meta, _ in variants],
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8", newline="\n")

    print(f"{len(variants)} simulations ({args.samples} samples, seed {args.seed}) -> {run_dir}")
    simulate(run_dir, args.samples, args.seed, args.jobs)
    analyze.analyze(run_dir)
    print(f"done: {run_dir / 'report.md'}")


if __name__ == "__main__":
    main()
