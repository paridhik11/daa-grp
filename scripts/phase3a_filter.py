#!/usr/bin/env python3
"""Phase 3A item 1: seed-selection cost variants and filter_auto versus filter / filter_fa.

Tuned and evaluated on E. coli 536, the synthetic genome and the stress text only (never on B. fragilis).
  A. selection variants (pure filter, 100 reads): original DP, prefix tables (all lengths), tables + length cap 1/2/3
  B. filter_auto vs filter vs filter_fa (fallback active), best-of-3 mean per setting
Writes results/phase3a_filter.md and results/phase3a_filter.csv.   Usage: phase3a_filter.py [--quick]
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, RES, run_bin  # noqa: E402
from simulate import make_reads, make_reads_lowcomplex  # noqa: E402

QUICK = "--quick" in sys.argv
NR = 20 if QUICK else 100
REP = 1 if QUICK else 3
ecoli_path = os.path.join(DATA, "ecoli536.txt")
ecoli = open(ecoli_path).read().strip()
md = ["## Phase 3A, item 1: cheaper seed selection and filter_auto\n"]
rows = []


def timed(algo, path, reads, k, env=None, nofb=False):
    old = {}
    for key, v in (env or {}).items():
        old[key] = os.environ.get(key); os.environ[key] = str(v)
    try:
        best = None
        for _ in range(REP):
            info, _s, tq = run_bin(algo, path, reads, k, nofallback=nofb, dump=False, tag="p3", times=True)
            ms = np.array([q["ns"] for q in tq]) / 1e6
            r = dict(mean=float(ms.mean()), med=float(np.median(ms)), sel_us=info["sel_ms"] / len(reads) * 1000, hits=info["hits"] / len(reads),
                     skip=info.get("auto_skip", 0), eq=info.get("auto_eq", 0), fa=info.get("auto_fa", 0), scan=info.get("auto_scan", 0))
            if best is None or r["mean"] < best["mean"]:
                best = r
        return best
    finally:
        for key, v in old.items():
            if v is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = v


# ---------------------------------------------------------------- A. selection variants
VARS = [("original DP (all lengths 1..11, per-seed lookups)", "filter_fa_full", {}),
        ("prefix tables, all lengths 1..11", "filter_fa", {"APPROX_SEL": 1}),
        ("tables + length cap 3", "filter_fa", {"APPROX_SEL": 2, "APPROX_CAP": 3}),
        ("tables + length cap 2 (default)", "filter_fa", {"APPROX_SEL": 2, "APPROX_CAP": 2}),
        ("tables + length cap 1", "filter_fa", {"APPROX_SEL": 2, "APPROX_CAP": 1})]
md.append("### A. Seed-selection variants (E. coli 536, m = 100, pure filter, mean over reads)\n")
md.append(f"Selection time is the mean time spent choosing seeds per query ({NR} reads, best of {REP} runs); hits are the mean seed hits of the chosen seeds, "
          "so a variant that is cheaper but chooses worse seeds shows up as more hits.\n")
md.append("| k | " + " | ".join(v[0] for v in VARS) + " |")
md.append("|---|" + "---|" * len(VARS))
KS_A = [1, 4, 8, 10, 12, 15] if not QUICK else [4, 8, 12]
for k in KS_A:
    reads = [r[0] for r in make_reads(ecoli, NR, 100, k, seed=11)]
    cells = []
    for name, algo, env in VARS:
        r = timed(algo, ecoli_path, reads, k, env, nofb=True)
        rows.append(dict(part="A", k=k, variant=name, sel_us=r["sel_us"], hits=r["hits"], mean_ms=r["mean"]))
        cells.append(f"{r['sel_us']:.1f} us, {r['hits']:.0f} hits, {r['mean']:.3f} ms")
    md.append(f"| {k} | " + " | ".join(cells) + " |")

# ---------------------------------------------------------------- B. filter_auto vs filter vs filter_fa
md.append("\n### B. filter_auto versus filter and filter_fa (fallback active; mean ms per query, best of %d runs)\n" % REP)
md.append("`auto/best` is the mean time of filter_auto divided by the smaller of the means of filter and filter_fa (target: at most 1.05). "
          "Choices are counts of queries per decision: skip = DP not run (a subset of eq), eq, fa, scan.\n")
md.append("| text | setting | filter | filter_fa | filter_auto | auto/best | choices (skip / eq / fa / scan) |")
md.append("|---|---|---|---|---|---|---|")
worst = []


def compare(label, path, reads, k, setting):
    a = timed("filter", path, reads, k)
    b = timed("filter_fa", path, reads, k)
    c = timed("filter_auto", path, reads, k)
    best = min(a["mean"], b["mean"])
    ratio = c["mean"] / best
    worst.append((ratio, label, setting))
    rows.append(dict(part="B", text=label, k=k, setting=setting, filter=a["mean"], filter_fa=b["mean"], auto=c["mean"], ratio=ratio))
    md.append(f"| {label} | {setting} | {a['mean']:.4f} | {b['mean']:.4f} | {c['mean']:.4f} | {ratio:.2f} | {c['skip']:.0f} / {c['eq']:.0f} / {c['fa']:.0f} / {c['scan']:.0f} |")


KS_B = [1, 2, 4, 6, 8, 10, 11, 12, 13, 14, 15, 17, 20] if not QUICK else [2, 8, 12, 17]
for k in KS_B:
    compare("E. coli 536", ecoli_path, [r[0] for r in make_reads(ecoli, NR, 100, k, seed=11)], k, f"m=100, k={k}")
for m, k in ([(200, 10), (400, 20), (800, 40)] if not QUICK else [(400, 20)]):
    compare("E. coli 536", ecoli_path, [r[0] for r in make_reads(ecoli, NR, m, k, seed=12)], k, f"m={m}, k={k}")
synth_path = os.path.join(DATA, "genome_synth.txt")
if os.path.exists(synth_path):
    synth = open(synth_path).read().strip()
    for k in ([2, 8, 10, 12, 15] if not QUICK else [10]):
        compare("synthetic", synth_path, [r[0] for r in make_reads(synth, NR, 100, k, seed=21)], k, f"m=100, k={k}")
sp, rp = os.path.join(DATA, "stress_text.txt"), os.path.join(DATA, "stress_regions.json")
if os.path.exists(sp) and os.path.exists(rp):
    st = open(sp).read().strip(); regions = [tuple(x) for x in json.load(open(rp))]
    for k in ((6, 10, 15) if not QUICK else (10,)):
        compare("stress", sp, [r[0] for r in make_reads_lowcomplex(st, regions, NR, 100, k, "edge", seed=34)], k, f"low-complexity edge, k={k}")
        compare("stress", sp, [r[0] for r in make_reads_lowcomplex(st, regions, NR, 100, k, "inside", seed=35)], k, f"low-complexity inside, k={k}")
        compare("stress", sp, [r[0] for r in make_reads(st, NR, 100, k, seed=32, mode="cluster")], k, f"clustered, k={k}")
worst.sort(reverse=True)
over = [w for w in worst if w[0] > 1.05]
md.append(f"\nWorst auto/best ratio: {worst[0][0]:.2f} ({worst[0][1]}, {worst[0][2]}); {len(over)} of {len(worst)} settings exceed 1.05"
          + ("" if not over else ": " + "; ".join(f"{w[1]} {w[2]} ({w[0]:.2f})" for w in over)) + ".")
open(os.path.join(RES, "phase3a_filter.md"), "w").write("\n".join(md) + "\n")
with open(os.path.join(RES, "phase3a_filter.csv"), "w", newline="") as f:
    keys = sorted({k for r in rows for k in r})
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
print("\n".join(md))
