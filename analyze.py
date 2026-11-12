"""Deterministic analysis of a sweep run folder.

Reads runs/<run>/manifest.json and runs/<run>/raw/*.csv (written by sweep.py
and qca_sim) and writes, into the same folder:

    segments.csv  one row per variant x segment
    summary.csv   one row per variant, with the verdict
    report.md     readable tables (removal, shift matrix, ...)

The same inputs always give byte-identical outputs: fixed rules, fixed
ordering, fixed number formatting.

    python analyze.py runs/2026-10-28_1730_all
"""
import argparse
import csv
import json
import sys
from pathlib import Path

INPUT_LABELS = ["A", "B", "C", "D", "F"]  # file labels; the paper calls F "E"
DIRECTIONS = ["up", "down", "left", "right"]

# Selvi et al. 2026, Table 4: single-cell omission, x = incorrect
PAPER_TABLE4_FAIL = {"01", "02", "04", "08"}


def read_raw(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["segment"] = int(r["segment"])
        r["bits"] = "".join(r[f"in_{k}"] for k in INPUT_LABELS)
        r["expected"] = int(sum(int(r[f"in_{k}"]) for k in INPUT_LABELS) >= 3)
        r["p"] = float(r["out_mid"])
        r["absmin"] = float(r["out_absmin"])
    return rows


def classify(manifest, run_dir):
    th = manifest["thresholds"]
    raw = run_dir / "raw"
    control = read_raw(raw / "control.csv")
    control_min = min(r["absmin"] for r in control)
    half = len(control) // 2

    seg_rows, summary = [], []
    for v in manifest["variants"]:
        rows = read_raw(raw / f"{v['name']}.csv")
        if len(rows) != len(control):
            sys.exit(f"{v['name']}: {len(rows)} segments, control has {len(control)}")
        wrong = []
        wrong_pass = [0, 0]
        for r, c in zip(rows, control):
            logic = int(r["p"] > 0)
            ok = logic == r["expected"] and abs(r["p"]) >= th["unreadable"]
            if not ok:
                wrong.append(r["bits"])
                wrong_pass[r["segment"] >= half] += 1
            seg_rows.append({
                "variant": v["name"], "segment": r["segment"], "pass": 1 + (r["segment"] >= half),
                "inputs_ABCDE": r["bits"], "expected": r["expected"], "out_p": f"{r['p']:.6f}",
                "out_absmin_hold": f"{r['absmin']:.6f}", "control_p": f"{c['p']:.6f}",
                "correct": int(ok), "same_as_control": int(logic == int(c["p"] > 0)),
            })
        min_abs = min(r["absmin"] for r in rows)
        if wrong:
            verdict = "FAIL"
        elif min_abs < th["weak_abs"] or min_abs < th["weak_rel"] * control_min:
            verdict = "WEAK"
        else:
            verdict = "OK"
        summary.append({
            "variant": v["name"], "type": v["type"],
            "paper_label": v.get("paper_label", ""), "file_index": v.get("file_index", ""),
            "role": v.get("role", ""), "direction": v.get("direction", ""),
            "distance_nm": fmt_num(v.get("distance_nm")),
            "min_gap_nm": fmt_num(v.get("min_gap_nm")),
            "overlap": "" if "overlap" not in v else int(v["overlap"]),
            "verdict": verdict, "mv5_correct": int(not wrong),
            "wrong_segments": len(wrong), "wrong_pass1": wrong_pass[0], "wrong_pass2": wrong_pass[1],
            "wrong_combos_ABCDE": " ".join(sorted(set(wrong))),
            "min_abs_p": f"{min_abs:.4f}", "control_min_abs_p": f"{control_min:.4f}",
            "rel_to_control": f"{min_abs / control_min:.3f}",
        })
    return seg_rows, summary


def fmt_num(x):
    return "" if x is None else f"{float(x):g}"


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def pct(n, d):
    return f"{100 * n / d:.2f} %" if d else "n/a"


def label_key(label):
    return (0, int(label)) if label.isdigit() else (1, label)


def build_report(manifest, summary):
    p = manifest["params"]
    L = [f"# Sweep report: {manifest['run_name']}", ""]
    L += ["## Run", "", md_table(["key", "value"], [
        ["kind", manifest["kind"]],
        ["created", manifest["created"]],
        ["design", f"`{manifest['design']}` (sha256 {manifest['hashes']['design'][:12]})"],
        ["qca_sim", f"sha256 {manifest['hashes']['qca_sim'][:12]}"],
        ["cells.csv", f"sha256 {manifest['hashes']['cells_csv'][:12]}"],
        ["engine", "bistable, exhaustive verification"],
        ["samples", p["samples"]],
        ["seed", p["seed"]],
        ["distances (nm)", " ".join(fmt_num(d) for d in p["distances"])],
        ["thresholds", ", ".join(f"{k}={v}" for k, v in manifest["thresholds"].items())],
        ["simulations", len(manifest["variants"])],
    ]), ""]

    ctrl = next(s for s in summary if s["type"] == "control")
    L += ["## Control", "",
          f"Verdict **{ctrl['verdict']}**, wrong segments {ctrl['wrong_segments']}, "
          f"min |P| in hold {ctrl['min_abs_p']}.", ""]

    removal = sorted((s for s in summary if s["type"] == "removal"), key=lambda s: label_key(s["paper_label"]))
    if removal:
        ok = sum(s["mv5_correct"] for s in removal)
        L += ["## Single-cell omission", "",
              f"Tolerated: **{ok}/{len(removal)} = {pct(ok, len(removal))}** "
              f"(paper Table 4: {len(removal) - len(PAPER_TABLE4_FAIL)}/{len(removal)} = "
              f"{pct(len(removal) - len(PAPER_TABLE4_FAIL), len(removal))}).", ""]
        rows = []
        for s in removal:
            ours = "MV5" if s["mv5_correct"] else "x"
            paper = "x" if s["paper_label"] in PAPER_TABLE4_FAIL else "MV5"
            rows.append([s["paper_label"], ours, paper, "yes" if ours == paper else "**no**",
                         f"{s['wrong_segments']}/{s['wrong_pass1']}/{s['wrong_pass2']}",
                         s["min_abs_p"], s["wrong_combos_ABCDE"] or "-"])
        L += [md_table(["cell", "ours", "paper", "match", "wrong (total/pass1/pass2)", "min abs P",
                        "wrong combos (ABCDE)"], rows), ""]

    shift = [s for s in summary if s["type"] == "shift"]
    if shift:
        dists = sorted({float(s["distance_nm"]) for s in shift})
        fails = [s for s in shift if not s["mv5_correct"]]
        L += ["## Single-cell displacement", "",
              f"Failing: **{len(fails)}/{len(shift)} = {pct(len(fails), len(shift))}**.", ""]

        rows = []
        for d in dists:
            ds = [s for s in shift if float(s["distance_nm"]) == d]
            df = [s for s in ds if not s["mv5_correct"]]
            gaps = sorted({float(s["min_gap_nm"]) for s in ds})
            rows.append([f"{d:g} nm", f"{len(df)}/{len(ds)}", pct(len(df), len(ds)),
                         sum(s["overlap"] == 1 for s in df), sum(s["overlap"] == 1 for s in ds),
                         " ".join(f"{g:g}" for g in gaps)])
        L += ["### By distance", "", md_table(
            ["distance", "fail", "fail %", "fail with overlap", "overlapping runs", "min gaps (nm)"], rows), ""]

        # Matrix in layout order (top-to-bottom, left-to-right)
        cells = sorted(manifest["cells"], key=lambda c: (c["y_nm"], c["x_nm"]))
        index = {(s["paper_label"], s["direction"], float(s["distance_nm"])): s for s in shift}
        header = ["cell", "role"] + [f"{d[0].upper()}{x:g}" for d in DIRECTIONS for x in dists] + ["fails"]
        rows = []
        for c in cells:
            row, n = [c["paper_label"], c["role"]], 0
            for d in DIRECTIONS:
                for x in dists:
                    s = index.get((c["paper_label"], d, x))
                    if s is None:
                        row.append("")
                        continue
                    n += not s["mv5_correct"]
                    row.append(f"{s['wrong_segments']}{'*' if s['overlap'] == 1 else ''}")
            rows.append(row + [f"{n}/{len(DIRECTIONS) * len(dists)}"])
        L += ["### Matrix", "",
              "Wrong segments out of 64 per cell, direction (U/D/L/R) and distance (nm). "
              "`*` = overlapping placement. Rows in layout order.", "",
              md_table(header, rows), ""]

        rows = [[s["variant"], s["paper_label"], s["direction"], s["distance_nm"], s["min_gap_nm"],
                 "yes" if s["overlap"] == 1 else "", s["wrong_segments"], s["min_abs_p"]] for s in fails]
        L += ["### Failing displacements", "", md_table(
            ["variant", "cell", "dir", "nm", "min gap", "overlap", "wrong", "min abs P"], rows), ""]

    faults = [s for s in summary if s["type"] != "control"]
    inconsistent = [s for s in faults if s["wrong_pass1"] != s["wrong_pass2"]]
    L += ["## Consistency between the two passes", "",
          "Every input combination is simulated twice (pass 1 = segments 0-31, pass 2 = 32-63).", ""]
    if inconsistent:
        L += [md_table(["variant", "wrong pass 1", "wrong pass 2"],
                       [[s["variant"], s["wrong_pass1"], s["wrong_pass2"]] for s in inconsistent]), ""]
    else:
        L += ["All variants give the same result in both passes.", ""]

    if faults:
        ps = [float(s["min_abs_p"]) for s in faults]
        weak = sum(s["verdict"] == "WEAK" for s in faults)
        L += ["## Output polarization", "",
              f"min |P| over all fault runs: {min(ps):.4f} (max {max(ps):.4f}); WEAK verdicts: {weak}.", ""]
    return "\n".join(L)


def analyze(run_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    seg_rows, summary = classify(manifest, run_dir)
    write_csv(run_dir / "segments.csv", seg_rows)
    write_csv(run_dir / "summary.csv", summary)
    report = build_report(manifest, summary)
    (run_dir / "report.md").write_text(report, encoding="utf-8", newline="\n")
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="folder with manifest.json and raw/")
    args = ap.parse_args()
    analyze(args.run_dir)
    print(f"wrote segments.csv, summary.csv, report.md in {args.run_dir}")


if __name__ == "__main__":
    main()
